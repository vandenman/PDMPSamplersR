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
  structure(list(
    n_observations = n_observations, subsample_size = subsample_size,
    prior_standata = prior_standata, residual_envelope = residual_envelope,
    anchor = if (is.null(anchor)) NULL else as.numeric(anchor)
  ), class = "stan_subsampling")
}

#' OMRF person-level residual envelope
#'
#' Constructs a stacked categorical-logit curvature bound for each person.
#' If `B_n` vertically stacks the node-specific affine logit maps, the block
#' diagonal categorical covariance has operator norm at most `1/2`, hence the
#' complete person Hessian is bounded by
#' `0.5 * opnorm(B_n)^2`. This is no larger than summing the separate node
#' spectral bounds and preserves a global, anchor-independent certificate.
#' Threshold and interaction arguments are Stan unconstrained block names and
#' are resolved at sampler initialization.
#'
#' @param X Integer person-by-node matrix encoded from zero.
#' @param seen Number of observed categories per node.
#' @param thresholds Threshold parameter block name.
#' @param interactions Interaction parameter block name; edge order is the
#'   upper triangle `(1,2), (1,3), ..., (P-1,P)`.
#' @return A certified Stan residual-envelope specification.
#' @rdname stan_subsampling
#' @export
omrf_residual_envelope <- function(X, seen, thresholds, interactions) {
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
      person_design[design_row + seq_len(q) - 1L, ] <- B
      design_row <- design_row + q
    }
    unique_weights[key_index] <- 0.5 * norm(person_design, type = "2")^2
  }
  weights <- unique_weights[match(row_keys, row_keys[unique_rows])]
  envelope <- stan_residual_envelope(weights)
  envelope$type <- "omrf"
  envelope$bound_type <- "stacked_person_spectral"
  envelope$X <- X
  envelope$seen <- seen
  envelope$thresholds <- thresholds
  envelope$interactions <- interactions
  envelope$edge_order <- edges
  class(envelope) <- c("omrf_residual_envelope", "stan_residual_envelope")
  envelope
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
#' @return A list of gradient-closure, residual-rate, envelope, and construction
#'   diagnostics.
#' @export
stan_subsampling_diagnostics <- function(path_to_stanmodel, standata,
                                    subsampling, position, subset,
                                    velocity = NULL,
                                    flow = c("ZigZag", "BouncyParticle",
                                             "AdaptiveBoomerang"),
                                    flow_mean = NULL, flow_cov = NULL) {
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
  result <- .pdmpsamplers_julia_eval(
    "PDMPSamplersRBridge.r_stan_subsampling_diagnostics(
      _diagnostic_model, _diagnostic_full_data, _diagnostic_prior_data,
      _diagnostic_subsampling, _diagnostic_position, _diagnostic_subset,
      _diagnostic_velocity, _diagnostic_flow, _diagnostic_flow_mean,
      _diagnostic_flow_cov
    );"
  )
  if (is.environment(result)) result <- as.list(result)
  result
}
