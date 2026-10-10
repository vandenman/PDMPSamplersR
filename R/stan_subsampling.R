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
      subsample_size < 1L || subsample_size > n_observations) {
    cli::cli_abort("Argument {.arg subsample_size} must satisfy 1 <= m <= N.")
  }
  n_observations <- as.integer(n_observations)
  subsample_size <- as.integer(subsample_size)
  if (missing(prior_standata) || is.null(prior_standata)) {
    cli::cli_abort("Argument {.arg prior_standata} is required.")
  }
  if (!inherits(residual_envelope, "stan_residual_envelope")) {
    cli::cli_abort("Argument {.arg residual_envelope} must be a certified package residual-envelope object.")
  }
  if (inherits(residual_envelope, "omrf_residual_envelope")) {
    n_persons <- residual_envelope$n_persons
    n_nodes <- residual_envelope$n_nodes
    factorization <- residual_envelope$factorization
    if (!rlang::is_integerish(n_persons, n = 1L, finite = TRUE) ||
        n_persons < 1L ||
        !rlang::is_integerish(n_nodes, n = 1L, finite = TRUE) ||
        n_nodes < 1L ||
        !is.character(factorization) || length(factorization) != 1L ||
        !factorization %in% c("person", "person_node")) {
      cli::cli_abort(paste0(
        "Specialized OMRF envelopes require positive integer {.field n_persons} ",
        "and {.field n_nodes} metadata and a supported {.field factorization}."))
    }
    expected <- if (factorization == "person") {
      as.integer(n_persons)
    } else {
      as.integer(n_persons) * as.integer(n_nodes)
    }
    if (n_observations != expected) {
      cli::cli_abort(paste0(
        "OMRF {.field ", factorization, "} factorization requires ", expected,
        " likelihood contributions; received ", n_observations, "."))
    }
    compact_norm <- factorization == "person_node" && isTRUE(
      residual_envelope$bound_type %in%
        c("node_local_norm_low_rank", "norm_local_hcv_compact"))
    if (compact_norm) {
      powers <- if (identical(residual_envelope$bound_type,
                              "norm_local_hcv_compact")) 4L else 3L
      weights <- residual_envelope$norm_weights
      if (!is.matrix(weights) || !is.numeric(weights) ||
          !identical(dim(weights), c(powers * as.integer(n_nodes),
                                     as.integer(n_persons))) ||
          any(!is.finite(weights)) || any(weights < 0)) {
        cli::cli_abort(paste0(
          "Compact OMRF norm weights must be a finite nonnegative ",
          powers, "P by N matrix matching {.field n_nodes} and ",
          "{.field n_persons}."))
      }
    }
    if (factorization == "person_node" &&
        residual_envelope$bound_type %in%
          c("node_box_first_order", "cell_marked_node_box_first_order")) {
      membership <- residual_envelope$node_box_membership
      nodes <- residual_envelope$node_box_nodes
      minima <- residual_envelope$node_box_minima
      maxima <- residual_envelope$node_box_maxima
      if (!is.matrix(membership) ||
          !identical(dim(membership), c(as.integer(n_nodes),
                                        as.integer(n_persons))) ||
          !is.numeric(membership) || any(!is.finite(membership)) ||
          any(membership != floor(membership)) || !length(nodes) ||
          any(membership < 1L | membership > length(nodes)) ||
          !is.matrix(minima) || !is.matrix(maxima) ||
          !identical(dim(minima), c(length(nodes), as.integer(n_nodes))) ||
          !identical(dim(maxima), dim(minima)) ||
          any(!is.finite(minima)) || any(!is.finite(maxima)) ||
          any(minima > maxima)) {
        cli::cli_abort(paste0(
          "Compact OMRF node-box specifications require P-by-N integer ",
          "membership and matching finite component-by-P extrema."))
      }
    }
  } else if (is.null(residual_envelope$weights) ||
             ncol(residual_envelope$weights) != n_observations) {
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
  # The OMRF simulation data can intentionally carry both independent-prior
  # and node-scale fields so that the same generated data object can be routed
  # to either Stan model.  For the native full-gradient independent model,
  # `prior_cauchy_scale` is the model-defining field and must take precedence;
  # node-scale prior-only data are stripped by `make_prior_stan_data()` before
  # this helper is called for the dependent route.
  interaction <- if (!is.null(standata$prior_cauchy_scale)) {
    list(type = "cauchy", scale = standata$prior_cauchy_scale)
  } else if (!is.null(standata$slab_base_scale) &&
             !is.null(standata$prior_global_scale_sd) &&
             !is.null(standata$prior_node_scale_sd)) {
    # Half-normal global and node scales (omrf_node_shared.stan).
    list(
      type = "node_scale_gaussian",
      scale = standata$slab_base_scale,
      global_scale_sd = standata$prior_global_scale_sd,
      node_scale_sd = standata$prior_node_scale_sd,
      node_scale_prior = "halfnormal"
    )
  } else if (!is.null(standata$slab_base_scale) &&
             !is.null(standata$prior_node_scale_sd)) {
    # No global scale and lognormal(0, prior_node_scale_sd) node scales.
    list(
      type = "node_scale_gaussian",
      scale = standata$slab_base_scale,
      node_scale_sd = standata$prior_node_scale_sd,
      node_scale_prior = "lognormal"
    )
  } else if (!is.null(standata$prior_interaction_sd)) {
    list(type = "gaussian", scale = standata$prior_interaction_sd)
  } else {
    cli::cli_abort(paste0(
      "The requested analytic OMRF prior is not a supported independent ",
      "Cauchy/Gaussian or node-scale Gaussian prior."))
  }
  values <- c(standata[required], interaction$scale,
              interaction$global_scale_sd, interaction$node_scale_sd)
  if (any(!vapply(values, function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) && value > 0
  }, logical(1)))) {
    cli::cli_abort("Analytic OMRF prior parameters must be finite and positive.")
  }
  spec$analytic_prior <- list(
    type = interaction$type,
    threshold_alpha = standata$prior_threshold_alpha,
    threshold_beta = standata$prior_threshold_beta,
    interaction_scale = interaction$scale,
    global_scale = if (!is.null(interaction$global_scale_sd))
      "global_scale" else NULL,
    node_scale = if (identical(interaction$type, "node_scale_gaussian"))
      "node_scale" else NULL,
    global_scale_sd = interaction$global_scale_sd,
    node_scale_sd = interaction$node_scale_sd,
    node_scale_prior = interaction$node_scale_prior)
  spec
}

#' OMRF person-level residual envelope
#'
#' Specifies exact subsampling for an ordinal Markov random field (OMRF)
#' likelihood with one likelihood contribution per person. The residual of a
#' person's gradient is taken relative to an analytic, undamped Hessian control
#' variate at the anchor, and every grid cell gets a certified roof on the
#' summed residual rates (geometry `"direct_factor_roof_person_hcv"`).
#' Threshold and interaction arguments are Stan unconstrained block names and
#' are resolved at sampler initialization.
#'
#' Several arguments accept a single value. They state the configuration
#' explicitly and are checked against it: `use_hcv = TRUE`,
#' `hcv_damping = Inf`, `hcv_warmup = "hcv"`, `factorization = "person"`,
#' `factor_sampling = "size_biased"`, `prior_backend = "analytic"`,
#' `exact_selected_refinement = TRUE`, `categorical_residual_cap = TRUE`, and
#' `geometry = "direct_factor_roof_person_hcv"`.
#'
#' @param X Integer person-by-node matrix encoded from zero.
#' @param seen Number of observed categories per node.
#' @param thresholds Threshold parameter block name.
#' @param interactions Interaction parameter block name; edge order is the
#'   upper triangle `(1,2), (1,3), ..., (P-1,P)`.
#' @param use_hcv,hcv_damping,hcv_warmup,factorization,factor_sampling,prior_backend,exact_selected_refinement,categorical_residual_cap,geometry
#'   Fixed configuration; see Details.
#' @param person_signed_roof Linear BPS certificate policy. Signed node
#'   residual endpoints are summed within each person before taking an
#'   absolute bound; the factor-summed certificate remains a valid fallback.
#'   The estimator and retained target are unchanged. Other flows keep the
#'   ordinary roof.
#' @param signed_only_person_roof Linear BPS policy which builds the certified
#'   signed whole-person endpoint roof without also building the factor-summed
#'   fallback table. It preserves the estimator and uses pointwise factor
#'   refinement for the selected person. Requires `person_signed_roof = TRUE`.
#' @param joint_signed_group_size Linear BPS joint signed marked roof. Use `4`
#'   for fixed seeded groups of four, `16`, `32` or `64` for weighted ragged
#'   groups, or `0` to keep the ordinary scheduler. Requires
#'   `signed_only_person_roof = TRUE`.
#' @param shared_current_categorical Evaluator that reuses current categorical
#'   moments for exact active-neighbour patterns. Requires joint groups.
#' @return A certified Stan residual-envelope specification.
#' @rdname stan_subsampling
#' @export
omrf_residual_envelope <- function(X, seen, thresholds, interactions,
                                   use_hcv = TRUE, hcv_damping = Inf,
                                   hcv_warmup = "hcv",
                                   factorization = "person",
                                   factor_sampling = "size_biased",
                                   prior_backend = "analytic",
                                   exact_selected_refinement = TRUE,
                                   categorical_residual_cap = TRUE,
                                   person_signed_roof = FALSE,
                                   signed_only_person_roof = FALSE,
                                   joint_signed_group_size = 0L,
                                   shared_current_categorical = FALSE,
                                   geometry = "direct_factor_roof_person_hcv") {
  fixed <- list(
    use_hcv = TRUE, hcv_damping = Inf, hcv_warmup = "hcv",
    factorization = "person", factor_sampling = "size_biased",
    prior_backend = "analytic", exact_selected_refinement = TRUE,
    categorical_residual_cap = TRUE,
    geometry = "direct_factor_roof_person_hcv")
  supplied <- list(
    use_hcv = use_hcv, hcv_damping = hcv_damping, hcv_warmup = hcv_warmup,
    factorization = factorization, factor_sampling = factor_sampling,
    prior_backend = prior_backend,
    exact_selected_refinement = exact_selected_refinement,
    categorical_residual_cap = categorical_residual_cap, geometry = geometry)
  for (name in names(fixed)) {
    if (!identical(supplied[[name]], fixed[[name]])) {
      cli::cli_abort("Argument {.arg {name}} only supports {.code {deparse(fixed[[name]])}}.")
    }
  }
  for (name in c("person_signed_roof", "signed_only_person_roof",
                 "shared_current_categorical")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      cli::cli_abort("Argument {.arg {name}} must be TRUE or FALSE.")
    }
  }
  if (signed_only_person_roof && !person_signed_roof) {
    cli::cli_abort(
      "{.arg signed_only_person_roof} requires {.arg person_signed_roof = TRUE}.")
  }
  if (!rlang::is_integerish(joint_signed_group_size, n = 1L,
                           finite = TRUE) ||
      !joint_signed_group_size %in% c(0L, 4L, 16L, 32L, 64L)) {
    cli::cli_abort("Argument {.arg joint_signed_group_size} must be 0, 4, 16, 32, or 64.")
  }
  if (joint_signed_group_size == 4L &&
      (!signed_only_person_roof || nrow(X) %% 4L != 0L)) {
    cli::cli_abort(paste0("Joint groups of four require a ",
      "divisible signed-only person-HCV population."))
  }
  if (shared_current_categorical && joint_signed_group_size == 0L) {
    cli::cli_abort("Shared-current categorical evaluation requires joint groups.")
  }
  if (joint_signed_group_size >= 16L &&
      !(shared_current_categorical && signed_only_person_roof &&
        nrow(X) >= 2L * joint_signed_group_size && all(seen == 4L))) {
    cli::cli_abort(paste0("Grouped joint signed bounds require N >= 2m, four categories per node, ",
      "and shared-current signed-only person roofs."))
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
  edges <- if (P > 1L) {
    do.call(rbind, lapply(seq_len(P - 1L), function(j) {
      cbind(j, seq.int(j + 1L, P))
    }))
  } else {
    matrix(integer(), nrow = 0L, ncol = 2L)
  }

  # The fields, their order and their values are part of the envelope
  # identity hashed into anchor-bank specifications, so the constant fields
  # below are kept although the Julia bridge supports only these values.
  factor_weights <- matrix(1, nrow = P, ncol = nrow(X))
  envelope <- stan_residual_envelope(factor_weights)
  envelope$type <- "omrf"
  envelope$bound_type <- "direct_factor_roof_person_hcv"
  envelope$person_signed_roof <- person_signed_roof
  if (signed_only_person_roof)
    envelope$signed_only_person_roof <- TRUE
  if (joint_signed_group_size > 0L)
    envelope$joint_signed_group_size <- as.integer(joint_signed_group_size)
  if (shared_current_categorical)
    envelope$shared_current_categorical <- TRUE
  envelope$X <- X
  envelope$seen <- seen
  envelope$thresholds <- thresholds
  envelope$interactions <- interactions
  envelope$backend <- "analytic"
  envelope$factorization <- "person"
  envelope$factor_sampling <- "size_biased"
  envelope$prior_backend <- "analytic"
  envelope$n_persons <- nrow(X)
  envelope$n_nodes <- P
  envelope$edge_order <- edges
  envelope$node_weights <- unname(factor_weights)
  envelope$covariate_compact <- FALSE
  envelope$use_hcv <- TRUE
  envelope$hcv_after_warmup <- FALSE
  envelope$hcv_damping <- Inf
  envelope$hcv_scope <- "global"
  envelope$hcv_bank_learning <- "geometric"
  envelope$zigzag_sign_aware <- FALSE
  envelope$exact_selected_refinement <- TRUE
  envelope$categorical_residual_cap <- TRUE
  envelope$signed_bps <- FALSE
  envelope$signed_harmonic <- FALSE
  envelope$direct_cv_cell_scope <- "cell"
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
