#' Compile and configure subsampling custom-Stan models
#'
#' `pdmp_subsample_hpp_path()` returns the installed thread-local subset-hook
#' header. `compile_pdmp_stan_model()` compiles a Stan source with the required
#' `--allow-undefined` and `USER_HEADER` settings.
#'
#' @param path Path to a Stan source file.
#' @return The header path, or the compiled shared-library path.
#' @name stan_subsampling
#' @export
pdmp_subsample_hpp_path <- function() {
  path <- system.file("stan", "pdmp_subsample.hpp", package = "PDMPSamplersR")
  if (!nzchar(path)) {
    candidate <- file.path("inst", "stan", "pdmp_subsample.hpp")
    if (file.exists(candidate)) path <- normalizePath(candidate)
  }
  if (!nzchar(path) || !file.exists(path)) {
    cli::cli_abort("The installed PDMP Stan subset header could not be found.")
  }
  path
}

#' @rdname stan_subsampling
#' @export
compile_pdmp_stan_model <- function(path) {
  validate_type(path, type = "character", n = 1)
  if (!file.exists(path) || !grepl("\\.stan$", path)) {
    cli::cli_abort("Argument {.arg path} must identify an existing Stan source file.")
  }
  check_for_julia_setup()
  .pdmpsamplers_julia_call(
    "_compile_model_with_header",
    normalizePath(path), pdmp_subsample_hpp_path()
  )
}

#' Certified residual envelope for a subsampling Stan likelihood
#'
#' @param weights Nonnegative component-by-observation curvature weights.
#' @param growth_rates Optional nonnegative growth coefficient per component.
#' @return A declarative residual-envelope specification.
#' @rdname stan_subsampling
#' @export
stan_residual_envelope <- function(weights, growth_rates = NULL) {
  if (is.vector(weights)) weights <- matrix(as.numeric(weights), nrow = 1L)
  if (!is.matrix(weights) || !is.numeric(weights) || length(weights) == 0L ||
      any(!is.finite(weights)) || any(weights < 0)) {
    cli::cli_abort("Argument {.arg weights} must be a nonempty finite nonnegative numeric matrix.")
  }
  if (is.null(growth_rates)) growth_rates <- rep(0, nrow(weights))
  if (!is.numeric(growth_rates) || length(growth_rates) != nrow(weights) ||
      any(!is.finite(growth_rates)) || any(growth_rates < 0)) {
    cli::cli_abort("Argument {.arg growth_rates} must contain one finite nonnegative value per envelope component.")
  }
  structure(
    list(type = "separable", weights = unname(weights),
         growth_rates = as.numeric(growth_rates)),
    class = "stan_residual_envelope"
  )
}

#' Exact observation-subsampling specification for a custom Stan model
#'
#' @param n_observations Total number of subsampling likelihood contributions.
#' @param subsample_size Number drawn without replacement at each proposal.
#' @param prior_standata Prior-only data list or JSON path for the same model.
#' @param residual_envelope A certified object from [stan_residual_envelope()]
#'   or a recognized package provider such as [omrf_residual_envelope()].
#' @param anchor Optional unconstrained anchor.
#' @return A validated observation-subsampling specification.
#' @rdname stan_subsampling
#' @export
stan_subsampling <- function(n_observations, subsample_size,
                                    prior_standata, residual_envelope,
                                    anchor = NULL) {
  if (!rlang::is_integerish(n_observations, n = 1L, finite = TRUE) ||
      n_observations < 2L) {
    cli::cli_abort("Argument {.arg n_observations} must be an integer of at least 2.")
  }
  if (!rlang::is_integerish(subsample_size, n = 1L, finite = TRUE) ||
      subsample_size < 1L || subsample_size >= n_observations) {
    cli::cli_abort("Argument {.arg subsample_size} must satisfy 1 <= m < N.")
  }
  n_observations <- as.integer(n_observations)
  subsample_size <- as.integer(subsample_size)
  if (missing(prior_standata) || is.null(prior_standata)) {
    cli::cli_abort("Argument {.arg prior_standata} is required.")
  }
  if (!inherits(residual_envelope, "stan_residual_envelope")) {
    cli::cli_abort("Argument {.arg residual_envelope} must be a certified package residual-envelope object.")
  }
  if (ncol(residual_envelope$weights) != n_observations) {
    cli::cli_abort("Residual-envelope weights must have one column per observation.")
  }
  if (!is.null(anchor) && (!is.numeric(anchor) || any(!is.finite(anchor)))) {
    cli::cli_abort("Argument {.arg anchor} must be NULL or a finite numeric vector.")
  }
  if (inherits(residual_envelope, "omrf_residual_envelope") &&
      identical(residual_envelope$prior_backend, "analytic")) {
    residual_envelope <- .attach_omrf_analytic_prior(
      residual_envelope, prior_standata, "prior_standata")
  }
  structure(list(
    n_observations = n_observations, subsample_size = subsample_size,
    prior_standata = prior_standata, residual_envelope = residual_envelope,
    anchor = if (is.null(anchor)) NULL else as.numeric(anchor)
  ), class = "stan_subsampling")
}

.attach_omrf_analytic_prior <- function(spec, standata, argument = "standata") {
  if (!is.list(standata)) {
    cli::cli_abort("Analytic OMRF priors require {.arg {argument}} as a named list.")
  }
  required <- c("prior_threshold_alpha", "prior_threshold_beta")
  if (!all(required %in% names(standata))) {
    cli::cli_abort("Analytic OMRF priors require threshold alpha and beta data.")
  }
  interaction <- if (!is.null(standata$prior_cauchy_scale)) {
    list(type = "cauchy", scale = standata$prior_cauchy_scale)
  } else if (!is.null(standata$prior_interaction_sd)) {
    list(type = "gaussian", scale = standata$prior_interaction_sd)
  } else {
    cli::cli_abort(paste0(
      "The requested analytic OMRF prior is not an independent Cauchy or ",
      "Gaussian prior."))
  }
  values <- c(standata[required], interaction$scale)
  if (any(!vapply(values, function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) && value > 0
  }, logical(1)))) {
    cli::cli_abort("Analytic OMRF prior parameters must be finite and positive.")
  }
  spec$analytic_prior <- list(
    type = interaction$type,
    threshold_alpha = standata$prior_threshold_alpha,
    threshold_beta = standata$prior_threshold_beta,
    interaction_scale = interaction$scale)
  spec
}

#' OMRF person-level residual envelope
#'
#' Constructs a stacked categorical-logit curvature bound for each person.
#' If `B_n` vertically stacks the node-specific affine logit maps, the block
#' diagonal categorical covariance has operator norm at most `1/2`, hence the
#' complete person Hessian is bounded by
#' `0.5 * opnorm(B_n)^2`. This is no larger than summing the separate node
#' spectral bounds and preserves a global, anchor-independent certificate.
#' The analytic backend uses a tighter low-rank node-local certificate. For
#' each node it expands the product of threshold and interaction terms against
#' the Euclidean norm of that person's neighbour-category vector. This needs
#' only three envelope components per node while retaining covariate magnitude
#' and preventing thresholds or nonincident interactions from inflating the
#' bound.
#' Threshold and interaction arguments are Stan unconstrained block names and
#' are resolved at sampler initialization.
#'
#' @param X Integer person-by-node matrix encoded from zero.
#' @param seen Number of observed categories per node.
#' @param thresholds Threshold parameter block name.
#' @param interactions Interaction parameter block name; edge order is the
#'   upper triangle `(1,2), (1,3), ..., (P-1,P)`.
#' @param backend Residual-gradient backend. The default `"analytic"` evaluates
#'   the exact OMRF likelihood residual without invoking Stan autodiff.
#'   `"stan"` retains the selected-likelihood Stan path for validation.
#' @param use_hcv Whether to use the analytic categorical Hessian control
#'   variate. This requires `backend = "analytic"` and retains a certified
#'   first-order fallback through damping.
#' @param hcv_damping Positive damping scale. By default this is the number of
#'   OMRF likelihood coordinates.
#' @param hcv_warmup HCV policy during warmup. `"hcv"` preserves the existing
#'   behavior. Experimental `"first_order"` uses the certified first-order
#'   residual bound during warmup and activates HCV at the main-phase boundary.
#' @param factorization Experimental likelihood contribution layout. `"person"`
#'   retains one contribution per person. `"person_node"` exposes each
#'   person-node conditional as a separate sparse contribution.
#' @param factor_sampling Factor-subset design. `"size_biased"` uses the
#'   generic contribution design; `"node_stratified"` draws an equal number
#'   from every node while retaining exact size-biased proposal mixing.
#' @param prior_backend Deterministic prior-gradient backend. `"stan"` supports
#'   arbitrary custom priors. `"analytic"` recognizes the package OMRF
#'   independent Gaussian and Cauchy prior contracts from `prior_standata`.
#' @param geometry Certified residual-rate geometry. `"norm_local"` uses a
#'   compact covariate-norm decomposition. `"covariate_local"` retains
#'   empirical neighbour covariates in threshold, edge, and edge-pair terms
#'   while sharing trajectory scales across people. `"pattern_local"` groups
#'   identical node-neighbour patterns and retains their signed directional
#'   logit geometry, which is especially tighter for scalar BPS/Boomerang
#'   rates but is more expensive to construct.
#' @return A certified Stan residual-envelope specification.
#' @rdname stan_subsampling
#' @export
omrf_residual_envelope <- function(X, seen, thresholds, interactions,
                                   backend = c("analytic", "stan"),
                                   use_hcv = FALSE, hcv_damping = NULL,
                                   hcv_warmup = c("hcv", "first_order"),
                                   factorization = c("person", "person_node"),
                                   factor_sampling = c("size_biased", "node_stratified"),
                                   prior_backend = c("stan", "analytic"),
                                   geometry = c("norm_local", "covariate_local",
                                                "pattern_local")) {
  backend <- match.arg(backend)
  hcv_warmup <- match.arg(hcv_warmup)
  factorization <- match.arg(factorization)
  factor_sampling <- match.arg(factor_sampling)
  prior_backend <- match.arg(prior_backend)
  geometry <- match.arg(geometry)
  if (!is.logical(use_hcv) || length(use_hcv) != 1L || is.na(use_hcv)) {
    cli::cli_abort("Argument {.arg use_hcv} must be TRUE or FALSE.")
  }
  if (isTRUE(use_hcv) && backend != "analytic") {
    cli::cli_abort("OMRF Hessian control variates require {.arg backend = \"analytic\"}.")
  }
  if (factorization == "person_node" && backend != "analytic") {
    cli::cli_abort("Person-node OMRF factorization requires {.arg backend = \"analytic\"}.")
  }
  if (factorization == "person" && factor_sampling != "size_biased") {
    cli::cli_abort("Node-stratified sampling requires {.arg factorization = \"person_node\"}.")
  }
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) || any(X != floor(X))) {
    cli::cli_abort("Argument {.arg X} must be a finite integer-valued person-by-node matrix.")
  }
  X <- matrix(as.integer(X), nrow(X), ncol(X))
  if (!rlang::is_integerish(seen, n = ncol(X), finite = TRUE) || any(seen < 1L)) {
    cli::cli_abort("Argument {.arg seen} must contain one positive category count per node.")
  }
  seen <- as.integer(seen)
  for (j in seq_len(ncol(X))) {
    if (any(X[, j] < 0L | X[, j] >= seen[j])) {
      cli::cli_abort("Column {j} of {.arg X} must be encoded in 0:(seen[j]-1).")
    }
  }
  if (!is.character(thresholds) || length(thresholds) != 1L ||
      !is.character(interactions) || length(interactions) != 1L) {
    cli::cli_abort("Arguments {.arg thresholds} and {.arg interactions} must be scalar block names.")
  }
  P <- ncol(X)
  threshold_starts <- cumsum(c(0L, seen[-length(seen)] - 1L))
  no_thresholds <- sum(seen - 1L)
  edges <- if (P > 1L) {
    do.call(rbind, lapply(seq_len(P - 1L), function(j) {
      cbind(j, seq.int(j + 1L, P))
    }))
  } else {
    matrix(integer(), nrow = 0L, ncol = 2L)
  }
  d_model <- no_thresholds + nrow(edges)
  if (is.null(hcv_damping)) hcv_damping <- d_model
  if (!is.numeric(hcv_damping) || length(hcv_damping) != 1L ||
      !is.finite(hcv_damping) || hcv_damping <= 0) {
    cli::cli_abort("Argument {.arg hcv_damping} must be one finite positive number.")
  }
  incidence <- lapply(seq_len(P), function(j) {
    incident <- which(edges[, 1L] == j | edges[, 2L] == j)
    list(
      edges = incident,
      neighbours = ifelse(edges[incident, 1L] == j,
                          edges[incident, 2L], edges[incident, 1L])
    )
  })
  row_keys <- apply(X, 1L, paste, collapse = ",")
  unique_rows <- which(!duplicated(row_keys))
  unique_weights <- numeric(length(unique_rows))
  unique_node_weights <- matrix(0, nrow = P, ncol = length(unique_rows))
  unique_hcv_remainder_weights <- matrix(0, nrow = P, ncol = length(unique_rows))
  for (key_index in seq_along(unique_rows)) {
    n <- unique_rows[key_index]
    person_design <- matrix(0, no_thresholds, d_model)
    design_row <- 1L
    for (j in seq_len(P)) {
      q <- seen[j] - 1L
      if (q == 0L) next
      B <- matrix(0, q, d_model)
      for (u in seq_len(q)) {
        B[u, threshold_starts[j] + u] <- 1
        B[u, no_thresholds + incidence[[j]]$edges] <-
          u * X[n, incidence[[j]]$neighbours]
      }
      Bnorm <- norm(B, type = "2")
      unique_node_weights[j, key_index] <- 0.5 * Bnorm^2
      # For p = softmax(z), C(p) = diag(p) - p p'. The identities
      # ||C(p)|| <= 1/2 and
      # ||C(p)-C(q)|| <= 3 ||p-q|| <= 3/2 ||z-w||
      # give the integrated categorical Taylor remainder
      # ||B'(p(Bx)-p(Ba))-B'C(p(Ba))B(x-a)||
      # <= (3/4) ||B||^3 ||x-a||^2.
      unique_hcv_remainder_weights[j, key_index] <- 0.75 * Bnorm^3
      person_design[design_row + seq_len(q) - 1L, ] <- B
      design_row <- design_row + q
    }
    unique_weights[key_index] <- 0.5 * norm(person_design, type = "2")^2
  }
  weights <- unique_weights[match(row_keys, row_keys[unique_rows])]
  node_weights <- unique_node_weights[, match(row_keys, row_keys[unique_rows]),
                                      drop = FALSE]
  hcv_remainder_weights <- unique_hcv_remainder_weights[
    , match(row_keys, row_keys[unique_rows]), drop = FALSE]
  neighbour_norms <- matrix(0, nrow = P, ncol = nrow(X))
  for (node in seq_len(P)) {
    neighbours <- incidence[[node]]$neighbours
    if (length(neighbours)) {
      neighbour_norms[node, ] <- sqrt(rowSums(
        X[, neighbours, drop = FALSE]^2))
    }
  }
  norm_weights <- matrix(0, nrow = 3L * P,
                         ncol = if (factorization == "person_node") nrow(X) * P else nrow(X))
  if (factorization == "person_node") {
    for (node in seq_len(P)) {
      columns <- (node - 1L) * nrow(X) + seq_len(nrow(X))
      norm_weights[node, columns] <- 1
      norm_weights[P + node, columns] <- neighbour_norms[node, ]
      norm_weights[2L * P + node, columns] <- neighbour_norms[node, ]^2
    }
  } else {
    norm_weights[seq_len(P), ] <- 1
    norm_weights[P + seq_len(P), ] <- neighbour_norms
    norm_weights[2L * P + seq_len(P), ] <- neighbour_norms^2
  }
  if (geometry == "covariate_local" && factorization != "person") {
    cli::cli_abort("Covariate-local geometry currently requires person factorization.")
  }
  covariate_weights <- matrix(numeric(), nrow = 0L, ncol = nrow(X))
  if (factorization == "person") {
    for (node in seq_len(P)) {
      neighbours <- incidence[[node]]$neighbours
      node_covariates <- X[, neighbours, drop = FALSE]
      node_rows <- matrix(1, nrow = 1L, ncol = nrow(X))
      if (length(neighbours)) {
        node_rows <- rbind(node_rows, t(node_covariates), t(node_covariates))
        for (left_edge in seq_along(neighbours)) {
          for (right_edge in seq_along(neighbours)) {
            node_rows <- rbind(node_rows,
              node_covariates[, left_edge] * node_covariates[, right_edge])
          }
        }
      }
      covariate_weights <- rbind(covariate_weights, node_rows)
    }
  }
  # The structural envelope groups factors that have the same node-local
  # neighbour configuration.  Its Julia scale callback uses the exact range
  # of the displaced categorical logits, rather than collapsing the design to
  # one spectral norm.  A factor belongs to exactly one component.
  component_nodes <- integer()
  component_covariates <- matrix(integer(), nrow = 0L, ncol = P)
  component_membership <- vector("list", P)
  for (node in seq_len(P)) {
    neighbours <- incidence[[node]]$neighbours
    keys <- if (length(neighbours)) {
      apply(X[, neighbours, drop = FALSE], 1L, paste, collapse = ",")
    } else {
      rep.int("", nrow(X))
    }
    unique_keys <- unique(keys)
    membership <- match(keys, unique_keys)
    component_membership[[node]] <- membership
    for (pattern in seq_along(unique_keys)) {
      representative <- which(membership == pattern)[1L]
      component_nodes <- c(component_nodes, node)
      component_covariates <- rbind(component_covariates, X[representative, ])
    }
  }
  component_offsets <- cumsum(c(0L, vapply(component_membership,
                                            max, integer(1L))))
  n_components <- tail(component_offsets, 1L)
  structural_groups <- if (factorization == "person_node") {
    matrix(integer(nrow(X) * P), nrow = 1L)
  } else {
    matrix(integer(P * nrow(X)), nrow = P)
  }
  for (node in seq_len(P)) {
    rows <- component_offsets[node] + component_membership[[node]]
    columns <- if (factorization == "person_node") {
      (node - 1L) * nrow(X) + seq_len(nrow(X))
    } else {
      seq_len(nrow(X))
    }
    if (factorization == "person_node") {
      structural_groups[1L, columns] <- rows
    } else {
      structural_groups[node, columns] <- rows
    }
  }
  if (factorization == "person_node") {
    factor_weights <- matrix(0, nrow = P, ncol = nrow(X) * P)
    factor_hcv_weights <- matrix(0, nrow = P, ncol = nrow(X) * P)
    for (node in seq_len(P)) {
      columns <- (node - 1L) * nrow(X) + seq_len(nrow(X))
      factor_weights[node, columns] <- node_weights[node, ]
      factor_hcv_weights[node, columns] <- hcv_remainder_weights[node, ]
    }
    envelope <- stan_residual_envelope(factor_weights)
    node_weights <- factor_weights
    hcv_remainder_weights <- factor_hcv_weights
  } else {
    envelope <- stan_residual_envelope(weights)
  }
  envelope$type <- "omrf"
  envelope$bound_type <- if (backend == "analytic" && !use_hcv &&
                             geometry == "pattern_local") {
    "pattern_local_range"
  } else if (backend == "analytic" && !use_hcv &&
             geometry == "covariate_local") {
    "covariate_local_expansion"
  } else if (backend == "analytic" && !use_hcv) {
    "node_local_norm_low_rank"
  } else {
    "stacked_person_spectral"
  }
  envelope$X <- X
  envelope$seen <- seen
  envelope$thresholds <- thresholds
  envelope$interactions <- interactions
  envelope$backend <- backend
  envelope$factorization <- factorization
  envelope$factor_sampling <- factor_sampling
  envelope$prior_backend <- prior_backend
  envelope$n_persons <- nrow(X)
  envelope$n_nodes <- P
  envelope$edge_order <- edges
  envelope$node_weights <- unname(node_weights)
  envelope$norm_weights <- unname(norm_weights)
  envelope$covariate_weights <- unname(covariate_weights)
  envelope$structural_groups <- unname(structural_groups)
  envelope$n_structural_components <- n_components
  envelope$component_nodes <- component_nodes
  envelope$component_covariates <- unname(component_covariates)
  envelope$use_hcv <- isTRUE(use_hcv)
  envelope$hcv_after_warmup <- isTRUE(use_hcv) &&
    identical(hcv_warmup, "first_order")
  envelope$hcv_damping <- as.numeric(hcv_damping)
  envelope$hcv_remainder_weights <- unname(hcv_remainder_weights)
  class(envelope) <- c("omrf_residual_envelope", "stan_residual_envelope")
  envelope
}

#' Native full-gradient backend for an independent-prior OMRF
#'
#' Constructs a native Julia backend for the complete, non-factorized OMRF
#' negative gradient and Hessian-vector product. It is intended for full
#' Adaptive Boomerang, Boomerang, BPS, and Zig-Zag samplers; it does not change
#' the sampler into factorized dynamics. The Stan model is still used for
#' unconstrained parameter metadata and output conversion.
#'
#' The model data passed to [pdmp_sample_from_stanmodel()] must contain the
#' independent threshold-prior parameters and either `prior_cauchy_scale` or
#' `prior_interaction_sd`. Threshold and interaction blocks must cover the
#' complete unconstrained state.
#'
#' @inheritParams omrf_residual_envelope
#' @return An OMRF full-gradient backend specification for the
#'   `full_gradient` argument of [pdmp_sample_from_stanmodel()].
#' @export
omrf_full_gradient <- function(X, seen, thresholds, interactions) {
  spec <- omrf_residual_envelope(
    X, seen, thresholds, interactions,
    backend = "analytic", use_hcv = FALSE,
    factorization = "person", factor_sampling = "size_biased",
    prior_backend = "analytic")
  structure(list(spec = spec), class = "omrf_full_gradient")
}

.resolve_unconstrained_spec <- function(spec, unc_names, arg = "parameters") {
  if (is.null(spec)) return(seq_along(unc_names))
  if (rlang::is_integerish(spec)) {
    idx <- as.integer(spec)
    if (any(idx < 1L | idx > length(unc_names))) {
      cli::cli_abort("Integer {.arg {arg}} values must lie in 1:{length(unc_names)}.")
    }
  } else if (is.character(spec)) {
    idx <- unlist(lapply(spec, function(value) {
      exact <- which(unc_names == value)
      if (length(exact)) return(exact)
      prefix <- which(startsWith(unc_names, paste0(value, ".")) |
                      startsWith(unc_names, paste0(value, "[")))
      if (!length(prefix)) {
        cli::cli_abort("{.arg {arg}} value {.val {value}} did not match an unconstrained name or block prefix.")
      }
      prefix
    }), use.names = FALSE)
  } else {
    cli::cli_abort("Argument {.arg {arg}} must be NULL, character, or integer.")
  }
  if (anyDuplicated(idx)) cli::cli_abort("Resolved {.arg {arg}} coordinates contain duplicates.")
  idx
}

#' Inspect Stan unconstrained parameter mappings
#'
#' @param path_to_stanmodel Stan source or compiled-library path.
#' @param standata Named data list or JSON path.
#' @param parameters Optional unconstrained names, block prefixes, or indices.
#' @return A data frame with resolved 1-based indices and names.
#' @rdname stan_subsampling
#' @export
stan_parameter_mapping <- function(path_to_stanmodel, standata,
                                   parameters = NULL) {
  validate_type(path_to_stanmodel, type = "character", n = 1)
  data_path <- standata
  temporary <- FALSE
  if (is.list(standata)) {
    data_path <- tempfile(fileext = ".json")
    write_stan_json(standata, data_path)
    temporary <- TRUE
  }
  if (temporary) on.exit(unlink(data_path), add = TRUE)
  check_for_julia_setup()
  model_path <- normalizePath(path_to_stanmodel, mustWork = TRUE)
  if (grepl("\\.stan$", model_path)) {
    code <- paste(readLines(model_path, warn = FALSE), collapse = "\n")
    model_path <- if (grepl("pdmp_get_subsample_", code, fixed = TRUE)) {
      compile_pdmp_stan_model(model_path)
    } else {
      .pdmpsamplers_julia_call("_compile_model", model_path)
    }
  }
  unc_names <- .pdmpsamplers_julia_call(
    "r_get_param_unc_names", model_path, normalizePath(data_path, mustWork = TRUE)
  )
  idx <- .resolve_unconstrained_spec(parameters, unc_names)
  data.frame(index = idx, name = unc_names[idx], stringsAsFactors = FALSE)
}

#' Diagnose a subsampling custom-Stan gradient and envelope
#'
#' Evaluates the persistent full, prior-only, and selected Stan modes at one
#' position without starting a sampler. The selected likelihood gradient is
#' unscaled; `selected_scale` reports the `N / m` factor applied by Julia when
#' forming the subsampling gradient. The residual event-rate diagnostic uses the
#' same typed envelope and trajectory geometry as sampling.
#'
#' @param path_to_stanmodel Stan source or compiled-library path containing the
#'   package subset hook.
#' @param standata Full-likelihood Stan data list or JSON path.
#' @param subsampling A specification from [stan_subsampling()].
#' @param position Finite unconstrained parameter vector.
#' @param subset One-based observation indices, with exactly the configured
#'   subsample size and no duplicates.
#' @param velocity Finite diagnostic velocity. Defaults to a vector of ones.
#' @param flow Supported subsampling flow used for the rate diagnostic.
#' @param flow_mean Optional flow mean.
#' @param flow_cov Optional flow covariance.
#' @param phase Diagnostic sampling phase. For staged HCV specifications,
#'   `"main"` activates HCV and `"warmup"` uses the first-order fallback.
#' @return A list of gradient-closure, residual-rate, envelope, and construction
#'   diagnostics.
#' @export
stan_subsampling_diagnostics <- function(path_to_stanmodel, standata,
                                    subsampling, position, subset,
                                    velocity = NULL,
                                    flow = c("ZigZag", "BouncyParticle",
                                             "AdaptiveBoomerang"),
                                    flow_mean = NULL, flow_cov = NULL,
                                    phase = c("main", "warmup")) {
  validate_type(path_to_stanmodel, type = "character", n = 1)
  if (!inherits(subsampling, "stan_subsampling")) {
    cli::cli_abort("Argument {.arg subsampling} must be created by {.fn stan_subsampling}.")
  }
  if (!is.numeric(position) || !length(position) || any(!is.finite(position))) {
    cli::cli_abort("Argument {.arg position} must be a nonempty finite numeric vector.")
  }
  if (!rlang::is_integerish(subset, n = subsampling$subsample_size,
                            finite = TRUE) ||
      any(subset < 1L | subset > subsampling$n_observations) ||
      anyDuplicated(subset)) {
    cli::cli_abort("Argument {.arg subset} must contain exactly m distinct indices in 1:N.")
  }
  if (is.null(velocity)) velocity <- rep(1, length(position))
  if (!is.numeric(velocity) || length(velocity) != length(position) ||
      any(!is.finite(velocity))) {
    cli::cli_abort("Argument {.arg velocity} must be finite and match {.arg position}.")
  }
  flow <- match.arg(flow)
  phase <- match.arg(phase)
  if (is.null(flow_mean)) flow_mean <- rep(0, length(position))
  if (is.null(flow_cov)) flow_cov <- diag(length(position))

  write_data <- function(data) {
    if (!is.list(data)) {
      validate_type(data, type = "character", n = 1)
      return(normalizePath(data, mustWork = TRUE))
    }
    path <- tempfile(fileext = ".json")
    write_stan_json(data, path)
    path
  }
  full_data <- write_data(standata)
  prior_data <- write_data(subsampling$prior_standata)
  temporary_data <- c(
    if (is.list(standata)) full_data else character(),
    if (is.list(subsampling$prior_standata)) prior_data else character()
  )
  if (length(temporary_data)) on.exit(unlink(temporary_data), add = TRUE)

  check_for_julia_setup()
  model_path <- normalizePath(path_to_stanmodel, mustWork = TRUE)
  if (grepl("\\.stan$", model_path)) {
    model_path <- compile_pdmp_stan_model(model_path)
  }
  JuliaCall::julia_assign("_diagnostic_model", model_path)
  JuliaCall::julia_assign("_diagnostic_full_data", full_data)
  JuliaCall::julia_assign("_diagnostic_prior_data", prior_data)
  JuliaCall::julia_assign("_diagnostic_subsampling", subsampling)
  JuliaCall::julia_assign("_diagnostic_position", as.numeric(position))
  JuliaCall::julia_assign("_diagnostic_subset", as.integer(subset))
  JuliaCall::julia_assign("_diagnostic_velocity", as.numeric(velocity))
  JuliaCall::julia_assign("_diagnostic_flow", flow)
  JuliaCall::julia_assign("_diagnostic_flow_mean", as.numeric(flow_mean))
  JuliaCall::julia_assign("_diagnostic_flow_cov", flow_cov)
  JuliaCall::julia_assign("_diagnostic_phase", phase)
  result <- .pdmpsamplers_julia_eval(
    "PDMPSamplersRBridge.r_stan_subsampling_diagnostics(
      _diagnostic_model, _diagnostic_full_data, _diagnostic_prior_data,
      _diagnostic_subsampling, _diagnostic_position, _diagnostic_subset,
      _diagnostic_velocity, _diagnostic_flow, _diagnostic_flow_mean,
      _diagnostic_flow_cov, _diagnostic_phase
    );"
  )
  if (is.environment(result)) result <- as.list(result)
  result
}
