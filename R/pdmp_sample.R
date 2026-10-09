with_seed <- function(code, seed = NULL) {
  if (is.null(seed)) {
    return(force(code))
  }

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, mode = "integer", inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, mode = "integer", inherits = FALSE) else NULL

  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, mode = "integer", inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(code)
}

.flow_spec_agreement <- function(resolved, requested_mean,
                                 requested_covariance, tolerance = 1e-12) {
  stopifnot(length(resolved) >= 1L)
  mean_agreement <- vapply(resolved, function(x)
    isTRUE(all.equal(as.numeric(x$mean), as.numeric(requested_mean),
                     tolerance = tolerance)), logical(1L))
  covariance_agreement <- vapply(resolved, function(x)
    isTRUE(all.equal(as.matrix(x$covariance), as.matrix(requested_covariance),
                     tolerance = tolerance, check.attributes = FALSE)),
    logical(1L))
  list(mean = mean_agreement, covariance = covariance_agreement,
       chain = mean_agreement & covariance_agreement)
}

.position_bank_result_spec <- function(counters, supplied, endpoint,
    selection_policy, unc_names, dataset_hash, envelope_hash, envelope,
    node_support_order) {
  counters <- if (is.list(counters)) counters else list()
  if (!is.null(supplied)) {
    # A supplied retained bank is immutable mathematical input. Counter
    # payloads are diagnostic and may be unavailable after a long JuliaCall
    # result conversion; they must not be the sole copy of restart state.
    if (!is.null(counters$anchor_positions) && !identical(
        unname(as.matrix(counters$anchor_positions)),
        unname(as.matrix(supplied$ordered_anchor_positions))))
      cli::cli_abort("The live position bank differs from the supplied frozen bank.")
    if (!is.null(counters$anchor_distance_scales) && !identical(
        as.numeric(counters$anchor_distance_scales),
        as.numeric(supplied$distance_scales)))
      cli::cli_abort("The live position-bank distance scales differ from the supplied frozen scales.")
    if (!is.null(supplied$warmup_bank_learning_policy) &&
        !is.null(counters$warmup_bank_learning_policy) &&
        !identical(as.character(counters$warmup_bank_learning_policy),
                   as.character(supplied$warmup_bank_learning_policy)))
      cli::cli_abort("The live HCV bank-learning policy differs from the supplied frozen bank.")
    spec <- supplied
    spec$immutable_envelope_grouping_spec <- envelope
  } else {
    if (is.null(counters$anchor_positions) ||
        is.null(counters$anchor_distance_scales))
      cli::cli_abort(paste0(
        "A newly prepared position bank was not returned by the live ",
        "sampler counters."))
    spec <- list(
      version = 1L,
      ordered_anchor_positions = unname(as.matrix(counters$anchor_positions)),
      distance_scales = as.numeric(counters$anchor_distance_scales),
      selection_policy = "position",
      tie_policy = "lowest_index_outgoing",
      requested_capacity = as.integer(counters$anchor_bank_requested_capacity),
      parameter_order = as.character(unc_names),
      dataset_identity = list(hash = dataset_hash),
      envelope_identity = list(hash = envelope_hash),
      immutable_envelope_grouping_spec = envelope)
    spec$warmup_bank_learning_policy <- as.character(
      counters$warmup_bank_learning_policy %||% "geometric")
    spec$warmup_bank_candidate_count <- as.integer(
      counters$warmup_bank_candidate_count %||% 0L)
    spec$warmup_bank_training_state_count <- as.integer(
      counters$warmup_bank_training_state_count %||% 0L)
    spec$warmup_bank_objective_single <- as.numeric(
      counters$warmup_bank_objective_single %||% NA_real_)
    spec$warmup_bank_objective_selected <- as.numeric(
      counters$warmup_bank_objective_selected %||% NA_real_)
    spec$warmup_bank_learning_seconds <- as.numeric(
      counters$warmup_bank_learning_seconds %||% 0)
  }
  spec$physical_endpoint <- endpoint
  spec$frozen_flow <- if (is.null(endpoint)) NULL else endpoint$flow
  if (identical(selection_policy, "node_position")) {
    spec$selection_policy <- "node_position"
    if (!is.null(counters$node_position_metric_variances))
      spec$physical_velocity_variances <-
        as.numeric(counters$node_position_metric_variances)
    spec$node_support_order <- node_support_order(envelope, unc_names)
    if (!is.null(counters$node_position_sources))
      spec$selected_sources <- as.integer(counters$node_position_sources)
    if (!is.null(counters$node_position_policy_version))
      spec$policy_version <- as.integer(counters$node_position_policy_version)
  }
  spec
}

# Internal function to validate and prepare common PDMP parameters
validate_pdmp_params <- function(d, flow, algorithm, T, t0 = 0.0, t_warmup = 0.0,
                                flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                                x0 = NULL, theta0 = NULL, show_progress = TRUE,
                                sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                                slab_prior = NULL,
                                grid_n = 30, grid_t_max = 2.0,
                                post_warmup_simplify = FALSE,
                                grid_bound = "constant",
                                grid_curvature_bound = NULL,
                                linear_area_threshold = 0.95,
                                linear_min_area_gain = 0.0,
                                n_chains = 1L, threaded = FALSE, seed = NULL,
                                adaptive_scheme = "diagonal",
                                unc_names = NULL,
                                lazy_low_tightness_threshold = 0.1,
                                lazy_max_low_tightness_rejections = 3L,
                                lazy_max_rejections = 0L,
                                allow_zero_adaptive_warmup = FALSE) {

  # Validate basic parameters
  d <- cast_integer(d, n = 1)
  validate_type(d, type = "integer", n = 1, positive = TRUE)
  validate_type(T, type = "double", n = 1, positive = TRUE)
  validate_type(t0, type = "double", n = 1)
  validate_type(t_warmup, type = "double", n = 1)

  if (t_warmup < 0)
    cli::cli_abort("Argument {.arg t_warmup} must be non-negative.")
  if (t_warmup >= T - t0)
    cli::cli_abort("Argument {.arg t_warmup} ({t_warmup}) must be less than {.code T - t0} ({T - t0}).")

  flow <- match.arg(flow, c("ZigZag", "BouncyParticle", "Boomerang", "AdaptiveBoomerang", "PreconditionedZigZag", "PreconditionedBPS"))
  grid_like_algorithms <- "GridThinningStrategy"
  algorithm <- match.arg(algorithm, c("ThinningStrategy", grid_like_algorithms, "RootsPoissonStrategy"))
  adaptive_scheme <- match.arg(adaptive_scheme, c("diagonal", "fullrank"))
  grid_bound <- match.arg(grid_bound, c("constant", "flat", "linear", "auto", "value_quadratic", "shared_node"))

  # AdaptiveBoomerang requires GridThinningStrategy and a warmup period
  if (flow == "AdaptiveBoomerang") {
    if (!algorithm %in% grid_like_algorithms)
      cli::cli_abort("{.val AdaptiveBoomerang} requires a grid-like algorithm.")
    if (t_warmup == 0 && !isTRUE(allow_zero_adaptive_warmup)) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for {.val AdaptiveBoomerang}.")
    }
  }

  # PreconditionedZigZag and PreconditionedBPS require GridThinningStrategy and a warmup period
  if (flow %in% c("PreconditionedZigZag", "PreconditionedBPS")) {
    if (!algorithm %in% grid_like_algorithms)
      cli::cli_abort("{.val {flow}} requires a grid-like algorithm.")
    if (t_warmup == 0 && !isTRUE(allow_zero_adaptive_warmup)) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for {.val {flow}}.")
    }
  }

  validate_type(c0, type = "double", n = 1, positive = TRUE)
  validate_type(show_progress, type = "logical", n = 1)

  n_chains <- cast_integer(n_chains, n = 1)
  validate_type(n_chains, type = "integer", n = 1, positive = TRUE)
  validate_type(threaded, type = "logical", n = 1)

  if (threaded && rlang::is_false(.pdmpsamplers_julia_eval("PDMPSamplersRBridge.r_threading_available()"))) {
    cli::cli_warn("Argument {.arg threaded} is set to TRUE but Julia was started with only one thread so this has no effect. Call Sys.setenv(\"JULIA_NUM_THREADS\"=<number>) to set the desired number of threads at the start of an analysis.")
  }

  if (!is.null(seed)) {
    if (!rlang::is_integerish(seed, n = 1)) {
      cli::cli_abort("Argument {.arg seed} must be NULL or an integerish scalar.")
    }
    seed <- as.integer(seed)
    if (seed < 0) {
      cli::cli_abort("Argument {.arg seed} must be non-negative.")
    }
  }

  # Handle and validate flow_mean and flow_cov
  if (is.null(flow_mean)) {
    flow_mean <- rep(0, d)
  } else {
    validate_type(flow_mean, type = "double", n = d)
  }

  if (is.null(flow_cov)) {
    flow_cov <- diag(1, nrow = d, ncol = d)
  } else {
    validate_type(flow_cov, type = "double", dims = c(d, d))
    if (!isSymmetric(flow_cov)) {
      cli::cli_abort("Argument {.arg flow_cov} must be a symmetric matrix.")
    }
  }

  # Validate x0
  if (is.null(x0)) {
    x0 <- with_seed(stats::rnorm(d), seed = seed)
  } else {
    validate_type(x0, type = "double", n = d)
  }

  # Validate theta0
  if (!is.null(theta0)) {
    validate_type(theta0, type = "double", n = d)
    if (flow == "ZigZag") {
      if (!all(theta0 %in% c(-1, 1))) {
        cli::cli_abort("Argument {.arg theta0} should only contain -1 or 1 when using ZigZag dynamics.")
      }
    }
  }

  # Validate sticky parameters
  validate_type(sticky, type = "logical", n = 1)
  if (!sticky && !is.null(slab_prior)) {
    cli::cli_abort("Argument {.arg slab_prior} requires {.arg sticky} to be {.code TRUE}.")
  }
  if (sticky) {
    if (!is.null(parameter_prior) && !is.null(slab_prior)) {
      cli::cli_abort("Use either legacy {.arg parameter_prior} or dependent {.arg slab_prior}, not both.")
    }
    if (!is.null(slab_prior) && !is.slab_prior(slab_prior)) {
      cli::cli_abort("Argument {.arg slab_prior} must be created by a dependent slab constructor.")
    }

    if (is.null(can_stick)) {
      if (!is.null(slab_prior)) {
        coef_idx <- .slab_coef_indices(slab_prior, d, unc_names = unc_names)
        if (is.null(coef_idx)) {
          cli::cli_abort("Dependent {.arg slab_prior} requires {.arg can_stick} unless {.arg coef} is supplied as integer parameter indices.")
        }
        can_stick <- rep(FALSE, d)
        can_stick[coef_idx] <- TRUE
      } else {
        can_stick <- rep(FALSE, d)
      }
    } else {
      validate_type(can_stick, type = "logical", n = d)
    }

    if (is.null(model_prior) || !is.model_prior(model_prior)) {
      cli::cli_abort("Argument {.arg model_prior} must be provided when {.arg sticky} is {.code TRUE}.")
    } else if (is.bernoulli(model_prior)) {
      if (is.null(slab_prior) && length(model_prior$prob) == 1) {
        model_prior <- bernoulli(prob = rep(model_prior$prob, d))
      } else if (is.null(slab_prior) && length(model_prior$prob) != d) {
        cli::cli_abort("For legacy sticky sampling, {.arg model_prior} of class {.cls bernoulli} must have a single probability or a vector of length {.arg d}.")
      } else if (!is.null(slab_prior)) {
        p <- .slab_beta_dimension(slab_prior, d, can_stick, unc_names = unc_names)
        if (!length(model_prior$prob) %in% c(1L, p, d)) {
          cli::cli_abort("For dependent {.arg slab_prior}, Bernoulli {.arg prob} must have length 1, the slab beta dimension ({p}), or {.arg d} ({d}).")
        }
      }
    }

    if (is.null(slab_prior)) {
      if (is.exchangeable_model_size_prior(model_prior)) {
        cli::cli_abort("{.fn exchangeable_model_size_prior} requires dependent {.arg slab_prior}; it is not supported by legacy sticky sampling.")
      }
      if (is.null(parameter_prior))
        cli::cli_abort("Argument {.arg parameter_prior} must be provided for legacy sticky sampling. Use {.arg slab_prior} for the dependent aggregate sticky path.")

      validate_type(parameter_prior, type = "double", n = d, positive = TRUE)
    } else {
      if (algorithm != "GridThinningStrategy") {
        cli::cli_abort("Dependent {.arg slab_prior} sticky sampling currently requires {.val GridThinningStrategy}.")
      }
      .validate_slab_prior_dimensions(slab_prior, d, can_stick, unc_names = unc_names, model_prior = model_prior)
    }

  }

  grid_n <- cast_integer(grid_n, n = 1)
  validate_type(grid_n, type = "integer", n = 1, positive = TRUE)
  validate_type(grid_t_max, type = "double", n = 1, positive = TRUE)
  validate_type(post_warmup_simplify, type = "logical", n = 1)
  if (!is.null(grid_curvature_bound)) {
    validate_type(grid_curvature_bound, type = "double", n = 1)
    if (!is.finite(grid_curvature_bound) || grid_curvature_bound < 0)
      cli::cli_abort("Argument {.arg grid_curvature_bound} must be NULL or a finite non-negative number.")
  }
  if (grid_bound == "shared_node" && is.null(grid_curvature_bound)) {
    cli::cli_abort("Argument {.arg grid_curvature_bound} is required for experimental {.val shared_node} bounds.")
  }
  validate_type(linear_area_threshold, type = "double", n = 1)
  validate_type(linear_min_area_gain, type = "double", n = 1)
  validate_type(lazy_low_tightness_threshold, type = "double", n = 1)
  lazy_max_low_tightness_rejections <- cast_integer(lazy_max_low_tightness_rejections, n = 1)
  lazy_max_rejections <- cast_integer(lazy_max_rejections, n = 1)
  validate_type(lazy_max_low_tightness_rejections, type = "integer", n = 1, positive = TRUE)
  validate_type(lazy_max_rejections, type = "integer", n = 1)
  if (linear_area_threshold < 0)
    cli::cli_abort("Argument {.arg linear_area_threshold} must be non-negative.")
  if (linear_min_area_gain < 0)
    cli::cli_abort("Argument {.arg linear_min_area_gain} must be non-negative.")
  if (lazy_low_tightness_threshold < 0)
    cli::cli_abort("Argument {.arg lazy_low_tightness_threshold} must be non-negative.")
  if (lazy_max_rejections < 0)
    cli::cli_abort("Argument {.arg lazy_max_rejections} must be non-negative.")

  return(list(
    d = d, flow = flow, algorithm = algorithm, T = T, t0 = t0, t_warmup = t_warmup,
    flow_mean = flow_mean, flow_cov = flow_cov, c0 = c0,
    x0 = x0, theta0 = theta0, show_progress = show_progress,
    sticky = sticky, can_stick = can_stick, model_prior = model_prior, parameter_prior = parameter_prior,
    slab_prior = slab_prior,
    grid_n = grid_n, grid_t_max = grid_t_max,
    post_warmup_simplify = post_warmup_simplify,
    grid_bound = grid_bound,
    grid_curvature_bound = grid_curvature_bound,
    linear_area_threshold = linear_area_threshold,
    linear_min_area_gain = linear_min_area_gain,
    lazy_low_tightness_threshold = lazy_low_tightness_threshold,
    lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
    lazy_max_rejections = lazy_max_rejections,
    n_chains = n_chains, threaded = threaded, seed = seed,
    adaptive_scheme = adaptive_scheme
  ))
}

.slab_coef_indices <- function(slab_prior, d, unc_names = NULL) {
  coef <- slab_prior$coef
  if (is.null(coef)) return(NULL)
  if (is.character(coef)) {
    if (is.null(unc_names)) {
      cli::cli_abort("Character slab {.arg coef} requires unconstrained parameter names; use integer indices for custom-gradient {.fn pdmp_sample}.")
    }
    return(.resolve_unconstrained_spec(coef, unc_names, "coef"))
  }
  idx <- as.integer(coef)
  if (any(idx < 1L) || any(idx > d)) {
    cli::cli_abort("Integer {.arg coef} values must be 1-based indices in 1:d.")
  }
  idx
}

.slab_beta_dimension <- function(slab_prior, d, can_stick, unc_names = NULL) {
  coef_idx <- .slab_coef_indices(slab_prior, d, unc_names = unc_names)
  if (is.null(coef_idx)) sum(can_stick) else length(coef_idx)
}

.slab_state_indices <- function(value, d, unc_names = NULL, arg = "logscale") {
  if (is.character(value)) {
    if (is.null(unc_names)) {
      cli::cli_abort("Character slab {.arg {arg}} requires unconstrained parameter names; use integer indices for custom-gradient {.fn pdmp_sample}.")
    }
    return(.resolve_unconstrained_spec(value, unc_names, arg))
  }
  idx <- as.integer(value)
  if (any(idx < 1L) || any(idx > d)) {
    cli::cli_abort("Integer {.arg {arg}} values must be 1-based indices in 1:d.")
  }
  idx
}

.validate_slab_prior_dimensions <- function(slab_prior, d, can_stick, unc_names = NULL, model_prior = NULL) {
  if (is.null(slab_prior)) return(invisible(TRUE))
  .validate_slab_prior_sampling_contract(slab_prior)
  coef_idx <- .slab_coef_indices(slab_prior, d, unc_names = unc_names)
  p <- if (is.null(coef_idx)) sum(can_stick) else length(coef_idx)
  if (p == 0L) {
    cli::cli_abort("Dependent {.arg slab_prior} requires at least one stickable coefficient.")
  }
  if (!is.null(coef_idx) && !all(can_stick[coef_idx])) {
    cli::cli_abort("Explicit slab {.arg coef} must be a subset of coordinates enabled by {.arg can_stick}.")
  }
  if (!is.null(model_prior) && is.exchangeable_model_size_prior(model_prior) && length(model_prior$omega) != p + 1L) {
    cli::cli_abort("For {.fn exchangeable_model_size_prior}, {.arg omega} must have length slab beta dimension + 1 ({p + 1L}).")
  }
  if (slab_prior$type == "dense_gaussian") {
    if (length(slab_prior$mean) != p) {
      cli::cli_abort("Length of {.arg mean} must match the number of slab coefficients.")
    }
    if (!all(dim(slab_prior$cov) == c(p, p))) {
      cli::cli_abort("Dimensions of {.arg cov} must match the number of slab coefficients.")
    }
  } else if (slab_prior$type == "independent_slab_density") {
    if (!(length(slab_prior$kappa) %in% c(1L, p))) {
      cli::cli_abort("Argument {.arg kappa} must be length 1 or match the number of slab coefficients.")
    }
  } else if (slab_prior$type == "exchangeable_gaussian") {
    if (slab_prior$u + p * slab_prior$v <= 0) {
      cli::cli_abort("Exchangeable slab covariance requires {.code u + p * v > 0}.")
    }
  } else if (slab_prior$type == "independent_logscale_gaussian") {
    logscale_idx <- .slab_state_indices(slab_prior$logscale, d, unc_names, arg = "logscale")
    if (length(logscale_idx) == 1L) logscale_idx <- rep(logscale_idx, p)
    if (length(logscale_idx) != p) {
      cli::cli_abort("Argument {.arg logscale} must have length 1 or match the number of slab coefficients.")
    }
    beta_idx <- if (is.null(coef_idx)) which(can_stick) else coef_idx
    if (length(slab_prior$log_base_scales) != 1L && length(slab_prior$log_base_scales) != p) {
      cli::cli_abort("Argument {.arg log_base_scales} must have length 1 or match the number of slab coefficients.")
    }
    if (length(intersect(beta_idx, logscale_idx)) > 0L) {
      cli::cli_abort("Slab {.arg logscale} coordinates must be disjoint from slab {.arg coef} coordinates.")
    }
    if (any(can_stick[unique(logscale_idx)])) {
      cli::cli_abort("Slab {.arg logscale} coordinates must be non-stickable.")
    }
  } else if (slab_prior$type == "loglinear_gaussian_scale") {
    logscale_idx <- .slab_state_indices(slab_prior$logscale, d, unc_names, arg = "logscale")
    beta_idx <- if (is.null(coef_idx)) which(can_stick) else coef_idx
    if (!identical(as.integer(slab_prior$logscale_design$dims),
                   c(as.integer(p), as.integer(length(logscale_idx))))) {
      cli::cli_abort("Argument {.arg logscale_design} must have one row per slab coefficient and one column per log-scale coordinate.")
    }
    if (!(length(slab_prior$log_base_scales) %in% c(1L, p))) {
      cli::cli_abort("Argument {.arg log_base_sd} must have length 1 or match the number of slab coefficients.")
    }
    if (anyDuplicated(logscale_idx)) {
      cli::cli_abort("Resolved {.arg logscale} coordinates cannot contain duplicates.")
    }
    if (length(intersect(beta_idx, logscale_idx)) > 0L) {
      cli::cli_abort("Slab {.arg logscale} coordinates must be disjoint from slab {.arg coef} coordinates.")
    }
    if (any(can_stick[logscale_idx])) {
      cli::cli_abort("Slab {.arg logscale} coordinates must be non-stickable.")
    }
  }
  invisible(TRUE)
}

.validate_slab_prior_sampling_contract <- function(slab_prior) {
  if (is.null(slab_prior)) return(invisible(TRUE))
  if (identical(slab_prior$type, "callback_gaussian") &&
      !rlang::is_function(slab_prior$active_prior_neggrad)) {
    cli::cli_abort(
      "{.fn gaussian_scale_mixture_slab} requires {.arg active_prior_neggrad} when used for dependent-slab sampling."
    )
  }
  invisible(TRUE)
}

.validate_dependent_slab_early <- function(slab_prior, sticky, flow, algorithm,
                                           model_prior, parameter_prior) {
  if (is.null(slab_prior)) return(invisible(TRUE))
  .validate_slab_prior_sampling_contract(slab_prior)
  if (!isTRUE(sticky)) {
    cli::cli_abort("Argument {.arg slab_prior} requires {.arg sticky} to be {.code TRUE}.")
  }
  if (!is.slab_prior(slab_prior)) {
    cli::cli_abort("Argument {.arg slab_prior} must be created by a dependent slab constructor.")
  }
  if (!is.null(parameter_prior)) {
    cli::cli_abort("Use either legacy {.arg parameter_prior} or dependent {.arg slab_prior}, not both.")
  }
  if (is.null(model_prior) || !is.model_prior(model_prior)) {
    cli::cli_abort("Argument {.arg model_prior} must be provided when {.arg sticky} is {.code TRUE}.")
  }
  if (algorithm != "GridThinningStrategy") {
    cli::cli_abort("Dependent {.arg slab_prior} sticky sampling currently requires {.val GridThinningStrategy}.")
  }
  invisible(TRUE)
}

#' Support Boundary Control
#'
#' Build a control list for support-boundary handling in PDMP samplers.
#'
#' Both \code{"line_search"} and \code{"line_search_truncated_refresh"} probe
#' along the linear ray \eqn{x_0 + t v} and are therefore only valid for
#' BPS/ZigZag-family flows with linear dynamics. For non-linear flows (e.g.,
#' Boomerang) these modes fall back to \code{"error"} behavior.
#'
#' @param mode Character string, how to handle support-boundary violations where
#'   the target or gradient becomes undefined during forward trajectory probing.
#'   One of \code{"error"} (default, fail fast), \code{"line_search"}
#'   (localize the first invalid time via bisection, then error), or
#'   \code{"line_search_truncated_refresh"} (heuristic BPS-family recovery that
#'   first searches for ordinary events before the localized boundary, handles
#'   such an event if one occurs first, and otherwise refreshes velocity at a
#'   valid interior point).
#' @param max_bisection_steps Integer, maximum bisection iterations.
#' @param time_rtol Numeric, relative tolerance for bisection.
#' @param time_atol Numeric, absolute tolerance for bisection.
#' @param clip_fraction Numeric in \code{(0, 1]}, fraction of the last-valid time
#'   used as a safe interior point after localization. The truncated-refresh mode
#'   may apply an additional conservative cap before refreshing.
#' @param max_refresh_attempts Integer, maximum number of refreshed velocities to
#'   try if \code{mode = "line_search_truncated_refresh"} reaches the support
#'   boundary before an ordinary event.
#' @param refresh_probe_time Numeric, short forward probe time used to reject
#'   immediately invalid refreshed velocities. A value of zero disables the
#'   probe; otherwise the Julia fallback may use a slightly longer scale-aware
#'   probe.
#' @param min_safe_time Numeric, minimum time gap used when clipping away from the
#'   localized boundary.
#'
#' @return A list suitable for the \code{support_boundary} argument.
#'
#' @export
support_boundary_control <- function(mode = c("error", "line_search", "line_search_truncated_refresh"),
                                     max_bisection_steps = 60L,
                                     time_rtol = 1e-8,
                                     time_atol = 1e-10,
                                     clip_fraction = 1 - 1e-10,
                                     max_refresh_attempts = 20L,
                                     refresh_probe_time = 1e-4,
                                     min_safe_time = 1e-12) {
  mode <- match.arg(mode)
  max_bisection_steps <- cast_integer(max_bisection_steps, n = 1)
  max_refresh_attempts <- cast_integer(max_refresh_attempts, n = 1)

  validate_type(max_bisection_steps, type = "integer", n = 1)
  validate_type(time_rtol, type = "double", n = 1)
  validate_type(time_atol, type = "double", n = 1)
  validate_type(clip_fraction, type = "double", n = 1, positive = TRUE)
  validate_type(max_refresh_attempts, type = "integer", n = 1)
  validate_type(refresh_probe_time, type = "double", n = 1)
  validate_type(min_safe_time, type = "double", n = 1)

  if (max_bisection_steps < 0)
    cli::cli_abort("Argument {.arg max_bisection_steps} must be non-negative.")
  if (time_rtol < 0)
    cli::cli_abort("Argument {.arg time_rtol} must be non-negative.")
  if (time_atol < 0)
    cli::cli_abort("Argument {.arg time_atol} must be non-negative.")
  if (clip_fraction > 1)
    cli::cli_abort("Argument {.arg clip_fraction} must be in (0, 1].")
  if (max_refresh_attempts < 0)
    cli::cli_abort("Argument {.arg max_refresh_attempts} must be non-negative.")
  if (refresh_probe_time < 0)
    cli::cli_abort("Argument {.arg refresh_probe_time} must be non-negative.")
  if (min_safe_time < 0)
    cli::cli_abort("Argument {.arg min_safe_time} must be non-negative.")

  structure(
    list(
      mode = mode,
      max_bisection_steps = max_bisection_steps,
      time_rtol = time_rtol,
      time_atol = time_atol,
      clip_fraction = clip_fraction,
      max_refresh_attempts = max_refresh_attempts,
      refresh_probe_time = refresh_probe_time,
      min_safe_time = min_safe_time
    ),
    class = "pdmp_support_boundary_control"
  )
}

validate_support_boundary_control <- function(support_boundary) {
  if (is.null(support_boundary))
    return(support_boundary_control())
  if (!is.list(support_boundary))
    cli::cli_abort("Argument {.arg support_boundary} must be a list created by {.fn support_boundary_control}.")

  defaults <- support_boundary_control()
  unknown <- setdiff(names(support_boundary), names(defaults))
  if (length(unknown) > 0) {
    cli::cli_abort(c(
      "Unknown field{?s} in {.arg support_boundary}.",
      "x" = "Got {.field {unknown}}."
    ))
  }

  do.call(support_boundary_control, utils::modifyList(defaults, support_boundary))
}


#' PDMP Sampling
#'
#' Performs Piecewise Deterministic Markov Process (PDMP) sampling using the
#' PDMPSamplers.jl Julia package.
#'
#' @param f Function that computes the NEGATIVE gradient. Should take a numeric vector
#'   and return a numeric vector of the same length.
#' @param d Integer, dimension of the problem.
#' @param flow Character string specifying the flow type. One of "ZigZag",
#'   "BouncyParticle", "Boomerang", "AdaptiveBoomerang", "PreconditionedZigZag",
#'   or "PreconditionedBPS".
#'   The `"AdaptiveBoomerang"` flow learns its reference (mean and precision)
#'   during warmup and requires a grid-like algorithm. It is supported with
#'   dependent slabs when `GridThinningStrategy` is selected; provide positive
#'   `t_warmup`, or let the sampler choose its 20% default.
#'   The `"PreconditionedZigZag"` and `"PreconditionedBPS"` flows learn a
#'   diagonal preconditioner during warmup and also require
#'   a grid-like algorithm.
#' @param algorithm Character string specifying the algorithm. One of
#'   "ThinningStrategy", "GridThinningStrategy", or "RootsPoissonStrategy".
#' @param T Numeric, total sampling time (default: 50000).
#' @param t0 Numeric, initial time (default: 0.0).
#' @param t_warmup Numeric, warmup time (default: 0.0). Events during warmup are discarded.
#' @param flow_mean Numeric vector of length d, mean vector for the flow (default: zero vector).
#' @param flow_cov Numeric matrix of size d x d, covariance matrix for the flow (default: identity matrix).
#' @param c0 Numeric, bound parameter (default: 1e-2).
#' @param x0 Numeric vector of length d, initial position (default: random normal).
#' @param theta0 Numeric vector of length d, initial velocity (default: random based on the flow).
#' @param hessian Function that returns the negative Hessian matrix of size d x d (default: NULL).
#'   Only used with GridThinningStrategy to compute the Hessian-vector product.
#' @param use_fd_hvp Logical compatibility switch selecting finite-difference
#'   curvature when no exact Hessian-vector product is used.
#' @param sticky Logical, whether to use sticky sampling (default: FALSE).
#' @param can_stick Logical vector of length d, which coordinates can stick (default: all FALSE).
#' @param model_prior Prior distribution object for model selection. Legacy
#'   sticky sampling accepts \code{bernoulli()} or \code{betabernoulli()}.
#'   Dependent-slab sampling also accepts \code{exchangeable_model_size_prior()}.
#' @param parameter_prior Numeric vector of length d, prior parameters for sticky sampling (default: NULL).
#' @param slab_prior Optional dependent slab prior created by
#'   \code{dense_gaussian_slab()}, \code{exchangeable_gaussian_slab()},
#'   \code{independent_slab_density()}, \code{gaussian_scale_mixture_slab()},
#'   or \code{arbitrary_slab_boundary()}. Mutually exclusive with
#'   \code{parameter_prior}. The full target supplied through \code{f} must
#'   contain the complete coherent slab represented by \code{slab_prior}; the
#'   package replaces that slab by its active-face restriction internally.
#' @param grid_n Integer, number of grid points for GridThinningStrategy (default: 30).
#' @param grid_t_max Numeric, maximum time for grid in GridThinningStrategy (default: 2.0).
#' @param post_warmup_simplify Logical. If \code{TRUE}, the grid-thinning
#'   strategy may switch to a constant bound after warmup, reducing
#'   gradient calls during the main sampling phase.
#' @param grid_bound Character string selecting the GridThinning proposal
#'   envelope. \code{"constant"} preserves the historical behavior.
#'   \code{"flat"}, \code{"linear"}, and \code{"auto"} expose signed-rate
#'   affine envelope machinery when supported by the chosen flow and derivative
#'   provider. \code{"shared_node"} is an experimental numerical clock using
#'   consecutive rate nodes and secant-disagreement inflation; it is not a
#'   mathematically certified thinning envelope.
#' @param grid_curvature_bound Optional nonnegative curvature bound required by
#'   the experimental `"shared_node"` grid bound.
#' @param linear_area_threshold Numeric gate for \code{grid_bound = "linear"}
#'   or \code{"auto"}; affine cells are skipped when their relative area gain is
#'   too small.
#' @param linear_min_area_gain Numeric absolute area-gain gate for affine cells.
#' @param lazy_low_tightness_threshold Nonnegative threshold used by lazy grid
#'   refinement diagnostics.
#' @param lazy_max_low_tightness_rejections Nonnegative number of consecutive
#'   low-tightness rejections allowed before refinement.
#' @param lazy_max_rejections Nonnegative total lazy-rejection limit; zero uses
#'   the package default behavior.
#' @param show_progress Logical, whether to show progress bar (default: TRUE).
#' @param n_chains Integer, number of chains to run (default: 1).
#' @param threaded Logical, whether to run chains in parallel (default: FALSE).
#' @param seed NULL (default) or a non-negative integer seed. The seed controls
#'   Julia's sampler RNG and any generated default initial position when
#'   \code{x0 = NULL}.
#' @param adaptive_scheme Character string, adaptation scheme for AdaptiveBoomerang.
#'   One of "diagonal" (default, O(d) per update) or "fullrank" (O(d^3) per update,
#'   better for correlated targets). Ignored for other flow types.
#' @param materialize Logical. If \code{TRUE} (default), the chain skeleton is
#'   immediately extracted from Julia into R so that \code{saveRDS()} /
#'   \code{readRDS()} work without any extra steps. Set to \code{FALSE} to skip
#'   the extraction and keep only the live Julia reference.
#'   This can save time and memory if you don't need to save the result or if you plan to call \code{materialize()} manually later.
#' @param online_summaries \code{NULL} (default) or options from
#'   \code{online_summary_options()}. Requires \code{trace_storage = "streaming"};
#'   the summaries are read with \code{read_online_summaries()}.
#' @param support_boundary A list created by \code{support_boundary_control()} that
#'   controls support-boundary diagnostics and heuristic event/refresh recovery.
#'
#' @return A \code{pdmp_result} object. Use \code{mean}, \code{var},
#'   \code{quantile}, etc. for continuous-time estimators, or
#'   \code{discretize} to obtain a sample matrix.
#'
#' @export
pdmp_sample <- function(f, d,
                        flow = c("ZigZag", "BouncyParticle", "Boomerang", "AdaptiveBoomerang", "PreconditionedZigZag", "PreconditionedBPS"),
                        algorithm = c("ThinningStrategy", "GridThinningStrategy",
                                      "RootsPoissonStrategy"),
                        T = 50000, t0 = 0.0, t_warmup = 0.0,
                        flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                        x0 = NULL, theta0 = NULL,
                        hessian = NULL,
                        sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                        slab_prior = NULL,
                        grid_n = 30, grid_t_max = 2.0,
                        use_fd_hvp = FALSE,
                        post_warmup_simplify = FALSE,
                        grid_bound = c("constant", "flat", "linear", "auto", "value_quadratic", "shared_node"),
                        grid_curvature_bound = NULL,
                        linear_area_threshold = 0.95,
                        linear_min_area_gain = 0.0,
                        lazy_low_tightness_threshold = 0.1,
                        lazy_max_low_tightness_rejections = 3L,
                        lazy_max_rejections = 0L,
                        show_progress = TRUE,
                        n_chains = 1L, threaded = FALSE, seed = NULL,
                        adaptive_scheme = c("diagonal", "fullrank"),
                        refresh_rate = 1.0,
                        materialize = TRUE,
                        trace_storage = c("dense", "streaming"),
                        streaming_trace_directory = NULL,
                        streaming_trace_buffer_events = 4096L,
                        online_summaries = NULL,
                        support_boundary = support_boundary_control()) {

  # Validate function argument (fail fast before Julia setup)
  if (!rlang::is_function(f)) {
    cli::cli_abort(c(
      "{.arg f} must be a function.",
      "x" = "You've supplied an object of {.cls {class(f)}}."
    ))
  }

  # Use common validation function
  trace_storage <- match.arg(trace_storage)
  restore_online_env <- .set_online_summary_env(online_summaries, trace_storage)
  on.exit(restore_online_env(), add = TRUE)
  if (identical(trace_storage, "streaming")) {
    if (is.null(streaming_trace_directory) ||
        !is.character(streaming_trace_directory) ||
        length(streaming_trace_directory) != 1L ||
        is.na(streaming_trace_directory) ||
        !nzchar(streaming_trace_directory)) {
      cli::cli_abort(paste0(
        "Argument {.arg streaming_trace_directory} must be one nonempty ",
        "path when {.arg trace_storage = \"streaming\"}."))
    }
    if (!rlang::is_integerish(streaming_trace_buffer_events, n = 1L,
                              finite = TRUE) ||
        streaming_trace_buffer_events < 1L) {
      cli::cli_abort(paste0(
        "Argument {.arg streaming_trace_buffer_events} must be a positive ",
        "integer."))
    }
  }
  streaming_trace_buffer_events <- as.integer(streaming_trace_buffer_events)
  params <- validate_pdmp_params(d, flow, algorithm, T, t0, t_warmup, flow_mean, flow_cov,
                                c0, x0, theta0, show_progress,
                                sticky, can_stick, model_prior, parameter_prior, slab_prior,
                                grid_n, grid_t_max, post_warmup_simplify,
                                grid_bound, grid_curvature_bound,
                                linear_area_threshold, linear_min_area_gain,
                                n_chains, threaded, seed,
                                adaptive_scheme = adaptive_scheme)

  # Test the function with a sample input
  tryCatch({
    test_output <- f(params$x0)
    if (!is.numeric(test_output) || length(test_output) != params$d) {
      cli::cli_abort("{.arg f} must return a numeric vector of length {params$d}, but returned length {length(test_output)}.")
    }
  }, error = function(e) {
    cli::cli_abort(c(
      "{.arg f} failed on test input.",
      "x" = conditionMessage(e)
    ))
  })

  # Test the hessian with a sample input
  if (!is.null(hessian)) {
    if (!rlang::is_function(hessian)) {
      cli::cli_abort(c(
        "{.arg hessian} must be a function or NULL.",
        "x" = "You've supplied an object of {.cls {class(hessian)}}."
      ))
    }
    tryCatch({
      test_output <- hessian(params$x0)
      if (!is.numeric(test_output) || !is.matrix(test_output) || !all(dim(test_output) == c(params$d, params$d))) {
        cli::cli_abort("{.arg hessian} must return a numeric matrix of dimensions {params$d}x{params$d}.")
      }
    }, error = function(e) {
      cli::cli_abort(c(
        "{.arg hessian} failed on test input.",
        "x" = conditionMessage(e)
      ))
    })
  }

  check_for_julia_setup()

  support_boundary <- validate_support_boundary_control(support_boundary)

  # Pass arguments to Julia
  for (nm in names(params))
    JuliaCall::julia_assign(nm, params[[nm]])

  JuliaCall::julia_assign("f", f)
  JuliaCall::julia_command("grad!(out, x) = out .= f(x);")
  JuliaCall::julia_assign("hessian_f", hessian)
  JuliaCall::julia_assign("support_boundary_mode", support_boundary$mode)
  JuliaCall::julia_assign("support_boundary_max_bisection_steps", support_boundary$max_bisection_steps)
  JuliaCall::julia_assign("support_boundary_time_rtol", support_boundary$time_rtol)
  JuliaCall::julia_assign("support_boundary_time_atol", support_boundary$time_atol)
  JuliaCall::julia_assign("support_boundary_clip_fraction", support_boundary$clip_fraction)
  JuliaCall::julia_assign("support_boundary_max_refresh_attempts", support_boundary$max_refresh_attempts)
  JuliaCall::julia_assign("support_boundary_refresh_probe_time", support_boundary$refresh_probe_time)
  JuliaCall::julia_assign("support_boundary_min_safe_time", support_boundary$min_safe_time)
  JuliaCall::julia_assign("streaming_trace_directory",
    if (identical(trace_storage, "streaming"))
      normalizePath(streaming_trace_directory, mustWork = FALSE) else NULL)
  JuliaCall::julia_assign("streaming_trace_buffer_events",
                          streaming_trace_buffer_events)
  result <- .pdmpsamplers_julia_eval("PDMPSamplersRBridge.r_pdmp_custom(
    grad!, d, x0, flow, algorithm, flow_mean, flow_cov;
    c0 = c0, grid_n = grid_n, grid_t_max = grid_t_max,
    post_warmup_simplify = post_warmup_simplify,
    grid_bound = grid_bound,
    grid_curvature_bound = grid_curvature_bound,
    linear_area_threshold = linear_area_threshold,
    linear_min_area_gain = linear_min_area_gain,
    t0 = t0, T = T, t_warmup = t_warmup,
    hessian = hessian_f,
    sticky = sticky, can_stick = can_stick,
    model_prior = model_prior, parameter_prior = parameter_prior,
    slab_prior = slab_prior,
    show_progress = show_progress, n_chains = n_chains, threaded = threaded,
    seed = seed,
    adaptive_scheme = adaptive_scheme,
    support_boundary_mode = support_boundary_mode,
    support_boundary_max_bisection_steps = support_boundary_max_bisection_steps,
    support_boundary_time_rtol = support_boundary_time_rtol,
    support_boundary_time_atol = support_boundary_time_atol,
    support_boundary_clip_fraction = support_boundary_clip_fraction,
    support_boundary_max_refresh_attempts = support_boundary_max_refresh_attempts,
    support_boundary_refresh_probe_time = support_boundary_refresh_probe_time,
    support_boundary_min_safe_time = support_boundary_min_safe_time,
    streaming_trace_directory = streaming_trace_directory,
    streaming_trace_buffer_events = streaming_trace_buffer_events
  );")
  if (is.environment(result)) result <- as.list(result)
  result <- new_pdmp_result(
    chains   = result$chains,
    stats    = result$stats,
    d        = result$d,
    n_chains = result$n_chains
  )
  if (materialize) result <- .materialize(result)
  result

}

#' PDMP Sampling from Stan Model
#'
#' Performs Piecewise Deterministic Markov Process (PDMP) sampling from a Stan model
#' using the PDMPSamplers.jl Julia package with BridgeStan integration.
#'
#' The gradient and Hessian-vector product are automatically derived from the Stan model
#' via the BridgeStan extension of PDMPSamplers.jl. This means you do not need to provide
#' a separate gradient or Hessian function.
#'
#' @param path_to_stanmodel Character, path to a Stan model file (.stan) or a
#'   compiled Stan model (.so/.dll/.dylib). If a .stan file is provided,
#'   BridgeStan will compile it automatically.
#' @param standata Either a character path to the Stan data file (JSON format),
#'   or a named list that will be written to a temporary JSON file via
#'   [write_stan_json()].
#' @param prior_stanmodel Compatibility argument. The full-gradient
#'   dependent-slab path does not use it. Custom-Stan subsampling accepts only
#'   the same model as \code{path_to_stanmodel}; use \code{prior_standata} to
#'   disable likelihood contributions.
#' @param prior_standata Optional prior-only Stan data path or named list used
#'   by \code{subsampling}; it is not needed for full-gradient slabs.
#' @param curvature_backend Character string selecting directional curvature:
#'   `"exact"` or `"finite_difference"`. Custom-Stan subsampling has no
#'   deterministic exact HVP and therefore always requires finite differences.
#'   Finite differences use
#'   shifted exact gradients but approximate their directional derivative and
#'   therefore have step-size and truncation error. `NULL` preserves the
#'   legacy `use_fd_hvp` mapping.
#' @param subsampling Optional exact subsampling custom-Stan specification
#'   created by [stan_subsampling()]. The selected Stan branch must add
#'   the unscaled sum of the installed observations; Julia applies `N / m`.
#' @param full_gradient Optional non-factorized native full-gradient backend.
#'   Currently [omrf_full_gradient()] supports independent-prior OMRF models
#'   and provides an exact native Hessian-vector product for Boomerang-family
#'   and other full samplers.
#' @param theta0 Optional initial velocity forwarded to both full-gradient and
#'   custom-Stan subsampling initialization.
#' @param subsampling_warmup Warmup gradient backend for subsampling runs.
#'   `"subsampled"` preserves the existing behavior. `"full"` uses the exact
#'   native OMRF full gradient and HVP during warmup, then switches to the
#'   configured subsampling model for retained sampling.
#' @param subsampling_anchor_updates Number of prepared OMRF anchors added from
#'   continuous-time warmup means. With full-gradient subsampling warmup, zero
#'   prepares only the exact terminal warmup anchor.
#' @param subsampling_use_anchor_bank Retain prepared OMRF anchors and select
#'   the nearest entry at event boundaries. If false, updates replace the
#'   active anchor. Full-gradient subsampling warmup can enable the bank with
#'   zero scheduled updates because it prepares a terminal anchor.
#' @param subsampling_anchor_bank_capacity Maximum prepared OMRF anchors kept
#'   per chain.
#' @param subsampling_anchor_selection_policy Research anchor-bank policy.
#'   `"adaptive"` preserves the existing event-history/envelope policy;
#'   `"position"` selects one frozen finite bank source by weighted position;
#'   research-only `"node_position"` selects a first-order source separately
#'   for every node using the frozen weighted-logit Gram rule. Global
#'   `"position"` supports linear Zig--Zag/BPS and harmonic Boomerang flows;
#'   `"node_position"` remains available only for linear BPS flows.
#' @param subsampling_position_anchor_bank_spec Optional version-1 mathematical
#'   position-bank specification saved on a previous result. Restoration
#'   rebuilds prepared provider quantities from its ordered positions and
#'   frozen distance scales; it requires zero warmup and matching dataset,
#'   parameter-order, and envelope/grouping identities.
#' @param subsampling_main_anchor_refresh_distance Positive Euclidean distance
#'   from the nearest prepared anchor that triggers a new exact OMRF anchor
#'   during main sampling. The default `Inf` disables main-phase preparation.
#' @param subsampling_node_source_check_distance Nonnegative scaled Euclidean
#'   movement required before reconsidering node-local warmup-bank sources.
#'   Zero preserves checking at every event boundary. Skipped checks retain
#'   the previous certified node sources and do not change the target.
#' @param subsampling_require_no_main_full_data If `TRUE`, fail the run when
#'   retained sampling performs any full-gradient evaluation, prepares an
#'   anchor, or traverses every OMRF person-node conditional in native Julia.
#'   Setup and warmup work are allowed. This is intended as a hard benchmark
#'   invariant for frozen-control-variate and direct-clock experiments.
#' @param warmup_adaptation_interval Nonnegative time between dynamics
#'   adaptation attempts during warmup. `NULL` uses ten equally spaced warmup
#'   windows (`t_warmup / 10`); zero attempts adaptation at every event. The
#'   same interval is used by full-gradient and subsampling runs.
#' @param boomerang_variance_floor_fraction Nonnegative opt-in research
#'   safeguard for adaptive diagonal Boomerang flows. At every adaptation
#'   update, named `global_scale` and `node_scale.*` coordinates retain at
#'   least this fraction of their preadaptation variance. Zero (the default)
#'   disables the safeguard.
#' @param flow_specification Role of an explicitly supplied `flow_mean` or
#'   `flow_cov`. `"initial"` permits legitimate warmup adaptation and records
#'   every chain's retained flow. `"frozen"` requires every retained chain to
#'   match the request and fails closed otherwise. `"auto"` uses `"initial"`
#'   when `t_warmup > 0` and `"frozen"` otherwise.
#' @inheritParams pdmp_sample
#'
#' @return A \code{pdmp_result} object. Use \code{mean}, \code{var},
#'   \code{quantile}, etc. for continuous-time estimators, or
#'   \code{discretize} to obtain a sample matrix.
#'
#' @export
pdmp_sample_from_stanmodel <- function(path_to_stanmodel, standata,
                        prior_stanmodel = NULL, prior_standata = NULL,
                        subsampling = NULL,
                        full_gradient = NULL,
                        subsampling_warmup = c("subsampled", "full"),
                        subsampling_anchor_updates = 0L,
                        subsampling_use_anchor_bank = FALSE,
                        subsampling_anchor_bank_capacity = 8L,
                        subsampling_anchor_selection_policy = c("adaptive", "position", "node_position"),
                        subsampling_position_anchor_bank_spec = NULL,
                        subsampling_main_anchor_refresh_distance = Inf,
                        subsampling_node_source_check_distance = 0,
                        subsampling_require_no_main_full_data = FALSE,
                        subsampling_diagnostic_window = 0,
                        flow = c("ZigZag", "BouncyParticle", "Boomerang", "AdaptiveBoomerang", "PreconditionedZigZag", "PreconditionedBPS"),
                        algorithm = c("ThinningStrategy", "GridThinningStrategy",
                                      "RootsPoissonStrategy"),
                        T = 50000, t0 = 0.0, t_warmup = 0.0,
                        warmup_adaptation_interval = NULL,
                        boomerang_variance_floor_fraction = 0,
                        flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                        flow_specification = c("auto", "initial", "frozen"),
                        x0 = NULL, theta0 = NULL, free0 = NULL,
                        stored_velocity0 = NULL,
                        sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                        slab_prior = NULL,
                        grid_n = 30, grid_t_max = 2.0,
                        use_fd_hvp = FALSE,
                        curvature_backend = NULL,
                        post_warmup_simplify = FALSE,
                        grid_bound = c("constant", "flat", "linear", "auto", "value_quadratic", "shared_node"),
                        grid_curvature_bound = NULL,
                        linear_area_threshold = 0.95,
                        linear_min_area_gain = 0.0,
                        lazy_low_tightness_threshold = 0.1,
                        lazy_max_low_tightness_rejections = 3L,
                        lazy_max_rejections = 0L,
                        show_progress = TRUE,
                        n_chains = 1L, threaded = FALSE, seed = NULL,
                        adaptive_scheme = c("diagonal", "fullrank"),
                        refresh_rate = 1.0,
                        materialize = TRUE,
                        trace_storage = c("dense", "streaming"),
                        streaming_trace_directory = NULL,
                        streaming_trace_buffer_events = 4096L,
                        online_summaries = NULL,
                        support_boundary = support_boundary_control()) {

  bridge_profile_enabled <- identical(
    Sys.getenv("OMRF_PROFILE_EXCLUSIVE", "0"), "1")
  bridge_profile_started <- proc.time()[["elapsed"]]
  bridge_profile_last <- bridge_profile_started
  bridge_profile <- list()
  bridge_profile_mark <- function(name) {
    if (!bridge_profile_enabled) return(invisible(NULL))
    now <- proc.time()[["elapsed"]]
    bridge_profile[[name]] <<- now - bridge_profile_last
    bridge_profile_last <<- now
    invisible(NULL)
  }

  # Validate file paths on R side before setting up Julia
  validate_type(path_to_stanmodel, type = "character", n = 1)
  flow <- match.arg(flow)
  trace_storage <- match.arg(trace_storage)
  restore_online_env <- .set_online_summary_env(online_summaries, trace_storage)
  on.exit(restore_online_env(), add = TRUE)
  if (identical(trace_storage, "streaming")) {
    if (is.null(streaming_trace_directory) ||
        !is.character(streaming_trace_directory) ||
        length(streaming_trace_directory) != 1L ||
        is.na(streaming_trace_directory) ||
        !nzchar(streaming_trace_directory)) {
      cli::cli_abort(paste0(
        "Argument {.arg streaming_trace_directory} must be one nonempty ",
        "path when {.arg trace_storage = \"streaming\"}."))
    }
    if (!rlang::is_integerish(streaming_trace_buffer_events, n = 1L,
                              finite = TRUE) ||
        streaming_trace_buffer_events < 1L) {
      cli::cli_abort(paste0(
        "Argument {.arg streaming_trace_buffer_events} must be a positive ",
        "integer."))
    }
  }
  streaming_trace_buffer_events <- as.integer(streaming_trace_buffer_events)
  algorithm <- match.arg(algorithm)
  adaptive_scheme <- match.arg(adaptive_scheme)
  flow_specification <- match.arg(flow_specification)
  subsampling_warmup <- match.arg(subsampling_warmup)
  supplied_flow_specification <- !is.null(flow_mean) || !is.null(flow_cov)
  resolved_flow_specification <- if (!supplied_flow_specification) "generated" else
    if (flow_specification != "auto") flow_specification else
      if (t_warmup > 0) "initial" else "frozen"
  if (!is.null(warmup_adaptation_interval)) {
    if (!is.numeric(warmup_adaptation_interval) ||
        length(warmup_adaptation_interval) != 1L ||
        !is.finite(warmup_adaptation_interval) ||
        warmup_adaptation_interval < 0) {
      cli::cli_abort(paste0(
        "Argument {.arg warmup_adaptation_interval} must be NULL or a ",
        "finite nonnegative number."))
    }
    warmup_adaptation_interval <- as.numeric(warmup_adaptation_interval)
  }
  if (!is.numeric(boomerang_variance_floor_fraction) ||
      length(boomerang_variance_floor_fraction) != 1L ||
      !is.finite(boomerang_variance_floor_fraction) ||
      boomerang_variance_floor_fraction < 0) {
    cli::cli_abort(paste0(
      "Argument {.arg boomerang_variance_floor_fraction} must be a finite ",
      "nonnegative number."))
  }
  boomerang_variance_floor_fraction <-
    as.numeric(boomerang_variance_floor_fraction)
  if (!rlang::is_integerish(n_chains, n = 1L, finite = TRUE) || n_chains < 1L) {
    cli::cli_abort("Argument {.arg n_chains} must be a positive integer.")
  }
  n_chains <- as.integer(n_chains)
  if (!is.null(subsampling) &&
      !inherits(subsampling, "stan_subsampling")) {
    cli::cli_abort("Argument {.arg subsampling} must be created by {.fn stan_subsampling}.")
  }
  if (!is.null(full_gradient) &&
      !inherits(full_gradient, "omrf_full_gradient")) {
    cli::cli_abort("Argument {.arg full_gradient} must be created by {.fn omrf_full_gradient}.")
  }
  if (!is.null(full_gradient) && !is.null(subsampling)) {
    cli::cli_abort("Arguments {.arg full_gradient} and {.arg subsampling} are mutually exclusive.")
  }
  if (subsampling_warmup == "full") {
    if (is.null(subsampling)) {
      cli::cli_abort("Argument {.arg subsampling_warmup = \"full\"} requires {.arg subsampling}.")
    }
    if (!inherits(subsampling$residual_envelope, "omrf_residual_envelope") ||
        !identical(subsampling$residual_envelope$backend, "analytic") ||
        !identical(subsampling$residual_envelope$prior_backend, "analytic")) {
      cli::cli_abort(paste0(
        "Full-gradient subsampling warmup requires an analytic OMRF residual ",
        "envelope with {.arg prior_backend = \"analytic\"}."))
    }
    if (t_warmup <= 0) {
      cli::cli_abort("Full-gradient subsampling warmup requires positive {.arg t_warmup}.")
    }
  }
  if (!is.null(full_gradient)) {
    # Preserve an explicitly supplied prior-only specification. This is
    # necessary for dependent OMRF data objects that contain both independent
    # Cauchy and node-scale fields; inferring from the full data would select
    # the wrong analytic branch.
    if (is.null(full_gradient$spec$analytic_prior)) {
      full_gradient$spec <- .attach_omrf_analytic_prior(
        full_gradient$spec, standata, "standata")
    }
  }
  if (!rlang::is_integerish(subsampling_anchor_updates, n = 1L, finite = TRUE) ||
      subsampling_anchor_updates < 0L) {
    cli::cli_abort("Argument {.arg subsampling_anchor_updates} must be a nonnegative integer.")
  }
  if (!is.logical(subsampling_use_anchor_bank) ||
      length(subsampling_use_anchor_bank) != 1L || is.na(subsampling_use_anchor_bank)) {
    cli::cli_abort("Argument {.arg subsampling_use_anchor_bank} must be TRUE or FALSE.")
  }
  if (!rlang::is_integerish(subsampling_anchor_bank_capacity, n = 1L, finite = TRUE) ||
      subsampling_anchor_bank_capacity < 1L) {
    cli::cli_abort("Argument {.arg subsampling_anchor_bank_capacity} must be a positive integer.")
  }
  subsampling_anchor_updates <- as.integer(subsampling_anchor_updates)
  subsampling_anchor_bank_capacity <- as.integer(subsampling_anchor_bank_capacity)
  subsampling_anchor_selection_policy <- match.arg(
    subsampling_anchor_selection_policy)
  if (subsampling_anchor_selection_policy %in% c("position", "node_position") &&
      !isTRUE(subsampling_use_anchor_bank)) {
    cli::cli_abort(paste0(
      "Argument {.arg subsampling_anchor_selection_policy = \"position\"} ",
      "requires {.arg subsampling_use_anchor_bank = TRUE}."))
  }
  resolved_flow <- match.arg(flow)
  if (identical(subsampling_anchor_selection_policy, "position") &&
      !resolved_flow %in% c("ZigZag", "PreconditionedZigZag",
        "BouncyParticle", "PreconditionedBPS", "Boomerang",
        "AdaptiveBoomerang")) {
    cli::cli_abort(paste0(
      "Argument {.arg subsampling_anchor_selection_policy = \"position\"} ",
      "requires linear Zig--Zag/BPS or harmonic Boomerang dynamics."))
  }
  if (identical(subsampling_anchor_selection_policy, "node_position") &&
      !resolved_flow %in% c("BouncyParticle", "PreconditionedBPS")) {
    cli::cli_abort(paste0(
      "Argument {.arg subsampling_anchor_selection_policy = \"node_position\"} ",
      "is available only for linear BPS flows."))
  }
  if (!is.null(subsampling_position_anchor_bank_spec) &&
      !(subsampling_anchor_selection_policy %in% c("position", "node_position"))) {
    cli::cli_abort(paste0(
      "Argument {.arg subsampling_position_anchor_bank_spec} requires ",
      "{.arg subsampling_anchor_selection_policy = \"position\"}."))
  }
  if (!is.null(subsampling_position_anchor_bank_spec) && t_warmup > 0) {
    cli::cli_abort("A restored position anchor bank requires {.arg t_warmup = 0}.")
  }
  if (!is.numeric(subsampling_main_anchor_refresh_distance) ||
      length(subsampling_main_anchor_refresh_distance) != 1L ||
      is.na(subsampling_main_anchor_refresh_distance) ||
      subsampling_main_anchor_refresh_distance <= 0) {
    cli::cli_abort(paste0(
      "Argument {.arg subsampling_main_anchor_refresh_distance} must be ",
      "positive or Inf."))
  }
  subsampling_main_anchor_refresh_distance <-
    as.numeric(subsampling_main_anchor_refresh_distance)
  if (!is.numeric(subsampling_node_source_check_distance) ||
      length(subsampling_node_source_check_distance) != 1L ||
      !is.finite(subsampling_node_source_check_distance) ||
      subsampling_node_source_check_distance < 0) {
    cli::cli_abort(paste0(
      "Argument {.arg subsampling_node_source_check_distance} must be ",
      "finite and nonnegative."))
  }
  subsampling_node_source_check_distance <-
    as.numeric(subsampling_node_source_check_distance)
  if (!is.numeric(subsampling_diagnostic_window) ||
      length(subsampling_diagnostic_window) != 1L ||
      is.na(subsampling_diagnostic_window) ||
      !is.finite(subsampling_diagnostic_window) ||
      subsampling_diagnostic_window < 0) {
    cli::cli_abort(paste0(
      "Argument {.arg subsampling_diagnostic_window} must be a finite ",
      "nonnegative number."))
  }
  subsampling_diagnostic_window <- as.numeric(subsampling_diagnostic_window)
  if (!is.logical(subsampling_require_no_main_full_data) ||
      length(subsampling_require_no_main_full_data) != 1L ||
      is.na(subsampling_require_no_main_full_data)) {
    cli::cli_abort(paste0(
      "Argument {.arg subsampling_require_no_main_full_data} must be ",
      "TRUE or FALSE."))
  }
  if (isTRUE(subsampling_use_anchor_bank) && subsampling_anchor_updates == 0L &&
      subsampling_warmup != "full" &&
      is.null(subsampling_position_anchor_bank_spec)) {
    cli::cli_abort(paste0(
      "Argument {.arg subsampling_use_anchor_bank = TRUE} requires positive ",
      "{.arg subsampling_anchor_updates}, unless {.arg subsampling_warmup = \"full\"} ",
      "prepares the terminal warmup anchor."))
  }
  if (is.finite(subsampling_main_anchor_refresh_distance) &&
      !isTRUE(subsampling_use_anchor_bank)) {
    cli::cli_abort(paste0(
      "Finite {.arg subsampling_main_anchor_refresh_distance} requires ",
      "{.arg subsampling_use_anchor_bank = TRUE}."))
  }
  if (is.null(subsampling) &&
      (subsampling_anchor_updates > 0L || isTRUE(subsampling_use_anchor_bank))) {
    cli::cli_abort("Subsampling anchor controls require {.arg subsampling}.")
  }
  if (subsampling_anchor_updates > 0L &&
      (!inherits(subsampling$residual_envelope, "omrf_residual_envelope") ||
       !identical(subsampling$residual_envelope$backend, "analytic"))) {
    cli::cli_abort("Subsampling anchor updates currently require an analytic OMRF residual envelope.")
  }
  if (subsampling_anchor_updates > 0L && t_warmup <= 0) {
    cli::cli_abort("Subsampling anchor updates require positive {.arg t_warmup}.")
  }
  if (!is.null(subsampling)) {
    if (!algorithm %in% c("GridThinningStrategy", "ThinningStrategy")) {
      cli::cli_abort("Subsampling custom-Stan sampling requires {.val GridThinningStrategy} or {.val ThinningStrategy}.")
    }
    if (is.null(prior_standata)) {
      prior_standata <- subsampling$prior_standata
    }
  }

  validate_type(use_fd_hvp, type = "logical", n = 1)
  if (!is.null(subsampling)) {
    analytic_omrf_hvp <- inherits(subsampling$residual_envelope,
      "omrf_residual_envelope") &&
      identical(subsampling$residual_envelope$backend, "analytic") &&
      identical(subsampling$residual_envelope$analytic_prior$type,
                "node_scale_gaussian") &&
      is.list(slab_prior) && identical(slab_prior$type,
                                      "loglinear_gaussian_scale")
    if (!is.null(curvature_backend)) {
      curvature_backend <- match.arg(curvature_backend, c("exact", "finite_difference"))
      if (curvature_backend == "exact" && !analytic_omrf_hvp) {
        cli::cli_abort(paste0(
          "This custom-Stan subsampling specification has no deterministic exact HVP; ",
          "use {.val finite_difference} curvature."))
      }
    }
    if (analytic_omrf_hvp) {
      if (identical(curvature_backend, "finite_difference") || isTRUE(use_fd_hvp)) {
        cli::cli_inform(paste0(
          "Using the exact composed HVP for analytic dependent OMRF subsampling; ",
          "finite-difference curvature is not supported for this combination."))
      }
      curvature_backend <- "exact"
    } else {
      curvature_backend <- curvature_backend %||% "finite_difference"
    }
    use_fd_hvp <- curvature_backend == "finite_difference"
  } else if (is.null(curvature_backend)) {
    curvature_backend <- if (isTRUE(use_fd_hvp)) "finite_difference" else "exact"
  } else {
    curvature_backend <- match.arg(curvature_backend, c("exact", "finite_difference"))
    if (isTRUE(use_fd_hvp) && curvature_backend != "finite_difference") {
      cli::cli_abort("Argument {.arg use_fd_hvp = TRUE} conflicts with {.arg curvature_backend = {curvature_backend}}.")
    }
  }
  support_boundary <- validate_support_boundary_control(support_boundary)
  if (!is.null(subsampling) && support_boundary$mode != "error") {
    cli::cli_abort("Custom-Stan subsampling currently requires {.arg support_boundary} mode {.val error}.")
  }

  if (is.list(standata)) {
    standata_path <- tempfile(fileext = ".json")
    write_stan_json(standata, standata_path)
    on.exit(unlink(standata_path), add = TRUE)
  } else {
    validate_type(standata, type = "character", n = 1)
    standata_path <- standata
  }

  prior_standata_path <- NULL
  if (!is.null(prior_standata)) {
    if (is.list(prior_standata)) {
      prior_standata_path <- tempfile(fileext = ".json")
      write_stan_json(prior_standata, prior_standata_path)
      on.exit(unlink(prior_standata_path), add = TRUE)
    } else {
      validate_type(prior_standata, type = "character", n = 1)
      prior_standata_path <- prior_standata
    }
  }
  bridge_profile_mark("stan_json_and_input_seconds")
  .validate_dependent_slab_early(
    slab_prior, sticky, flow, algorithm, model_prior, parameter_prior
  )
  if (is.null(prior_stanmodel)) {
    prior_stanmodel <- path_to_stanmodel
  } else {
    validate_type(prior_stanmodel, type = "character", n = 1)
  }

  if (!is.null(subsampling) && is.null(prior_standata_path)) {
    cli::cli_abort("Subsampling custom-Stan sampling requires prior-only Stan data.")
  }

  if (!file.exists(path_to_stanmodel))
    cli::cli_abort("Stan model file not found: {.path {path_to_stanmodel}}")
  if (!file.exists(standata_path))
    cli::cli_abort("Stan data file not found: {.path {standata_path}}")
  if (!grepl("\\.(so|dll|dylib|stan)$", path_to_stanmodel))
    cli::cli_abort(c(
      "{.arg path_to_stanmodel} should point to a Stan model ({.file .stan}) or a compiled Stan model ({.file .so}, {.file .dll}, or {.file .dylib}).",
      "i" = "Got: {.path {path_to_stanmodel}}"
    ))
  if (!grepl("\\.json$", standata_path))
    cli::cli_abort(c(
      "{.arg standata} should be a JSON file path or a list that can be written to JSON.",
      "i" = "Got: {.path {standata_path}}",
      "i" = "Use {.fn write_stan_json} to create a data file or pass a list directly."
    ))
  if (!is.null(prior_standata_path)) {
    if (!file.exists(prior_stanmodel))
      cli::cli_abort("Prior Stan model file not found: {.path {prior_stanmodel}}")
    if (!file.exists(prior_standata_path))
      cli::cli_abort("Prior Stan data file not found: {.path {prior_standata_path}}")
    if (!grepl("\\.(so|dll|dylib|stan)$", prior_stanmodel))
      cli::cli_abort("{.arg prior_stanmodel} should point to a Stan model or compiled Stan library.")
    if (!grepl("\\.json$", prior_standata_path))
      cli::cli_abort("{.arg prior_standata} should be a JSON file path or a named list.")
  }

  if (!is.null(subsampling) &&
      !identical(normalizePath(prior_stanmodel, mustWork = TRUE),
                 normalizePath(path_to_stanmodel, mustWork = TRUE))) {
    cli::cli_abort(c(
      "Custom-Stan subsampling requires {.arg prior_stanmodel} to be the same model as {.arg path_to_stanmodel}.",
      "i" = "Use prior-only data to switch off likelihood contributions."
    ))
  }

  check_for_julia_setup()
  bridge_profile_mark("julia_setup_and_validation_seconds")

  # Normalize paths to absolute
  path_to_stanmodel <- normalizePath(path_to_stanmodel, mustWork = TRUE)
  standata_path     <- normalizePath(standata_path,     mustWork = TRUE)
  if (!is.null(prior_standata_path)) {
    prior_stanmodel <- normalizePath(prior_stanmodel, mustWork = TRUE)
    prior_standata_path <- normalizePath(prior_standata_path, mustWork = TRUE)
  }

  compile_control_enabled <- any(tolower(Sys.getenv(c(
    "PDMPSAMPLERSR_BRIDGESTAN_CACHE",
    "PDMPSAMPLERSR_BRIDGESTAN_NATIVE",
    "PDMPSAMPLERSR_BRIDGESTAN_STANC_O1"
  ), "false")) %in% c("1", "true", "yes", "y"))
  if (!is.null(subsampling) && grepl("\\.stan$", path_to_stanmodel)) {
    path_to_stanmodel <- .pdmpsamplers_julia_call(
      "_compile_model_with_header", path_to_stanmodel,
      pdmp_subsample_hpp_path()
    )
    prior_stanmodel <- path_to_stanmodel
  } else if (compile_control_enabled && grepl("\\.stan$", path_to_stanmodel)) {
    path_to_stanmodel <- .pdmpsamplers_julia_call("_compile_model", path_to_stanmodel)
  }
  if (is.null(subsampling) && compile_control_enabled && !is.null(prior_standata_path) && grepl("\\.stan$", prior_stanmodel)) {
    prior_stanmodel <- .pdmpsamplers_julia_call("_compile_model", prior_stanmodel)
  }
  bridge_profile_mark("stan_compile_seconds")

  JuliaCall::julia_assign("_path_to_stan_model", path_to_stanmodel)
  JuliaCall::julia_assign("_path_to_stan_data",  standata_path)
  JuliaCall::julia_assign("_path_to_prior_stan_data_full", prior_standata_path)
  JuliaCall::julia_assign("_subsampling", subsampling)
  JuliaCall::julia_assign("_full_gradient_spec",
                          if (is.null(full_gradient)) NULL else full_gradient$spec)
  JuliaCall::julia_assign("_n_chains", as.integer(n_chains))

  # Initialize exactly the persistent contexts that sampling will reuse, and
  # obtain dimension/name metadata from those contexts.
  JuliaCall::julia_assign("_curvature_backend", curvature_backend)
  if (is.null(subsampling)) {
    JuliaCall::julia_command("_stan_model = BridgeStan.StanModel(_path_to_stan_model, _path_to_stan_data; warn=false);")
    if (is.null(full_gradient)) {
      JuliaCall::julia_command("_pdmp_model = PDMPModel(_stan_model; hvp = _curvature_backend != \"finite_difference\");")
    } else {
      JuliaCall::julia_command("_pdmp_model = PDMPSamplersRBridge.build_omrf_full_model(_full_gradient_spec, BridgeStan.param_unc_names(_stan_model), Int(BridgeStan.param_unc_num(_stan_model)));")
    }
    d <- JuliaCall::julia_eval("_pdmp_model.d")
    unc_names <- JuliaCall::julia_eval("BridgeStan.param_unc_names(_stan_model)")
  } else {
    JuliaCall::julia_command("_prepared_stan_subsampling = PDMPSamplersRBridge.prepare_stan_subsampling(_path_to_stan_model, _path_to_stan_data, _path_to_prior_stan_data_full, _subsampling, _n_chains);")
    d <- JuliaCall::julia_eval("_prepared_stan_subsampling.d")
    unc_names <- JuliaCall::julia_eval("_prepared_stan_subsampling.unc_names")
  }
  bridge_profile_mark("julia_context_initialization_seconds")
  if (is.character(can_stick)) {
    stick_idx <- .resolve_unconstrained_spec(can_stick, unc_names, "can_stick")
    can_stick <- rep(FALSE, d)
    can_stick[stick_idx] <- TRUE
  }
  if (!is.null(subsampling)) {
    envelope_spec <- subsampling$residual_envelope
    if (inherits(envelope_spec, "omrf_residual_envelope")) {
      threshold_idx <- .resolve_unconstrained_spec(
        envelope_spec$thresholds, unc_names, "thresholds")
      interaction_idx <- .resolve_unconstrained_spec(
        envelope_spec$interactions, unc_names, "interactions")
      expected_thresholds <- sum(envelope_spec$seen - 1L)
      expected_interactions <- ncol(envelope_spec$X) *
        (ncol(envelope_spec$X) - 1L) / 2L
      if (length(threshold_idx) != expected_thresholds) {
        cli::cli_abort("Resolved OMRF threshold block has {length(threshold_idx)} coordinates; expected {expected_thresholds}.")
      }
      if (length(interaction_idx) != expected_interactions) {
        cli::cli_abort("Resolved OMRF interaction block has {length(interaction_idx)} coordinates; expected {expected_interactions}.")
      }
      if (length(intersect(threshold_idx, interaction_idx))) {
        cli::cli_abort("Resolved OMRF threshold and interaction blocks overlap.")
      }
    }
  }
  if (!is.null(full_gradient)) {
    envelope_spec <- full_gradient$spec
    threshold_idx <- .resolve_unconstrained_spec(
      envelope_spec$thresholds, unc_names, "thresholds")
    interaction_idx <- .resolve_unconstrained_spec(
      envelope_spec$interactions, unc_names, "interactions")
    prior_idx <- integer()
    if (identical(envelope_spec$analytic_prior$type,
                  "node_scale_gaussian")) {
      prior_idx <- c(
        if (!is.null(envelope_spec$analytic_prior$global_scale))
          .resolve_unconstrained_spec(
            envelope_spec$analytic_prior$global_scale, unc_names,
            "global_scale"),
        .resolve_unconstrained_spec(
          envelope_spec$analytic_prior$node_scale, unc_names,
          "node_scale"))
    }
    covered <- sort(c(threshold_idx, interaction_idx, prior_idx))
    if (!identical(covered, seq_len(d))) {
      cli::cli_abort(c(
        "Native OMRF full gradients require threshold and interaction blocks to cover the complete unconstrained state.",
        "i" = "Covered {length(covered)} of {d} coordinates."
      ))
    }
  }

  # Use common validation function
  params <- validate_pdmp_params(d, flow, algorithm, T, t0, t_warmup, flow_mean, flow_cov,
                                 c0, x0, theta0, show_progress,
                                 sticky, can_stick, model_prior, parameter_prior, slab_prior,
                                 grid_n, grid_t_max, post_warmup_simplify,
                                 grid_bound, grid_curvature_bound,
                                 linear_area_threshold, linear_min_area_gain,
                                 n_chains, threaded, seed,
                                 adaptive_scheme = adaptive_scheme,
                                 unc_names = unc_names,
                                 lazy_low_tightness_threshold = lazy_low_tightness_threshold,
                                 lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
                                 lazy_max_rejections = lazy_max_rejections,
                                 allow_zero_adaptive_warmup =
                                   !is.null(subsampling_position_anchor_bank_spec))
  if (!is.null(free0)) {
    if (!isTRUE(sticky) || length(free0) != d || anyNA(free0))
      cli::cli_abort("{.arg free0} must be a complete logical sticky-state mask of length {d}.")
    free0 <- as.logical(free0)
  }
  if (!is.null(stored_velocity0)) {
    if (!isTRUE(sticky) || length(stored_velocity0) != d ||
        any(!is.finite(stored_velocity0)))
      cli::cli_abort("{.arg stored_velocity0} must be a finite sticky-state vector of length {d}.")
    stored_velocity0 <- as.numeric(stored_velocity0)
  }
  if (xor(is.null(free0), is.null(stored_velocity0)))
    cli::cli_abort("{.arg free0} and {.arg stored_velocity0} must be supplied together for continuation.")

  if (!is.null(subsampling) && isTRUE(params$threaded) && params$n_chains > 1L) {
    cli::cli_warn(c(
      "Custom-Stan subsampling serializes every BridgeStan call across Julia threads to avoid Stan Math autodiff memory corruption.",
      "i" = "This prevents segmentation faults but may limit parallel-chain speedups.",
      "i" = "For true parallel speedups with Stan-backed models, use separate R/Julia processes or a Julia-native thread-safe gradient implementation."
    ))
  }

  verify_flow_specification <- supplied_flow_specification
  flow_hash <- .pdmp_object_hash
  # Flow construction may make an algebraically neutral covariance ->
  # precision -> scale -> covariance round trip.  Hash the numerical
  # representation at the same precision used by the fail-closed comparison,
  # rather than treating sub-ULP reconstruction noise as different provenance.
  flow_vector_hash <- function(value) flow_hash(signif(as.numeric(value), 12L))
  flow_matrix_hash <- function(value) {
    value <- as.matrix(value)
    flow_hash(list(dimensions = as.integer(dim(value)),
                   values = signif(as.numeric(value), 12L)))
  }
  dataset_identity_hash <- .pdmp_dataset_hash
  if (!is.null(subsampling_position_anchor_bank_spec)) {
    specs_to_check <- if (!is.null(subsampling_position_anchor_bank_spec$version))
      list(subsampling_position_anchor_bank_spec) else
      subsampling_position_anchor_bank_spec
    dataset_hash <- dataset_identity_hash(standata)
    # The categorical residual cap changes only retained thinning roofs; it
    # does not change bank positions, grouping, metric, or preparation.  Keep
    # it outside the mathematical bank identity so a frozen bank can be used
    # in the declared before/after certificate intervention.
    bank_envelope_spec <- subsampling$residual_envelope
    bank_envelope_spec$categorical_residual_cap <- NULL
    envelope_hash <- flow_hash(bank_envelope_spec)
    for (spec in specs_to_check) {
      if (!identical(as.integer(spec$version), 1L) ||
          !identical(spec$selection_policy, subsampling_anchor_selection_policy) ||
          !identical(spec$tie_policy, "lowest_index_outgoing")) {
        cli::cli_abort("Unsupported position anchor bank specification.")
      }
      if (identical(spec$selection_policy, "node_position") &&
          (!identical(as.integer(spec$policy_version), 1L) ||
           is.null(spec$physical_velocity_variances) ||
           is.null(spec$node_support_order))) {
        cli::cli_abort("Incomplete node-position mathematical bank specification.")
      }
      dataset_ok <- identical(spec$dataset_identity$hash, dataset_hash)
      envelope_ok <- identical(spec$envelope_identity$hash, envelope_hash)
      parameter_order_ok <- identical(as.character(spec$parameter_order),
                                      as.character(unc_names))
      if (!dataset_ok || !envelope_ok || !parameter_order_ok) {
        if (identical(Sys.getenv("OMRF_PROFILE_EXCLUSIVE", "0"), "1") &&
            !envelope_ok) {
          message("OMRF bank envelope diagnostic: expected=",
                  spec$envelope_identity$hash, " actual=", envelope_hash,
                  " comparison=", paste(all.equal(
                    spec$immutable_envelope_grouping_spec,
                    subsampling$residual_envelope), collapse = "; "))
        }
        cli::cli_abort(paste0(
          "The restored position anchor bank does not match the dataset, ",
          "parameter order, or immutable envelope/grouping specification ",
          "(dataset=", dataset_ok, ", parameter_order=", parameter_order_ok,
          ", envelope=", envelope_ok, ")."))
      }
      if (!is.null(spec$frozen_flow)) {
        same_bank_mean <- isTRUE(all.equal(
          as.numeric(spec$frozen_flow$mean), as.numeric(params$flow_mean),
          tolerance = 1e-12, check.attributes = FALSE))
        same_bank_covariance <- isTRUE(all.equal(
          as.matrix(spec$frozen_flow$covariance),
          as.matrix(params$flow_cov), tolerance = 1e-12,
          check.attributes = FALSE))
        if (!same_bank_mean || !same_bank_covariance) {
          cli::cli_abort(paste0(
            "The restored position anchor bank was prepared under a ",
            "different frozen flow mean or covariance."))
        }
      }
    }
  }
  requested_flow_provenance <- NULL
  if (verify_flow_specification) {
    resolved_preflight <- .pdmpsamplers_julia_call(
      "r_resolve_flow_spec", params$flow, params$flow_mean, params$flow_cov,
      adaptive_scheme = params$adaptive_scheme,
      refresh_rate = as.numeric(refresh_rate))
    if (is.environment(resolved_preflight)) resolved_preflight <- as.list(resolved_preflight)
    same_mean <- isTRUE(all.equal(as.numeric(resolved_preflight$mean),
                                  as.numeric(params$flow_mean), tolerance = 1e-12))
    same_cov <- isTRUE(all.equal(as.matrix(resolved_preflight$covariance),
                                 as.matrix(params$flow_cov), tolerance = 1e-12,
                                 check.attributes = FALSE))
    if (!same_mean || !same_cov) {
      cli::cli_abort(c(
        "The resolved Julia flow does not match the requested flow specification.",
        "i" = "Mean agreement: {same_mean}; covariance agreement: {same_cov}.",
        "i" = "Refusing to begin retained sampling."
      ))
    }
    requested_flow_provenance <- list(
      specification_role = resolved_flow_specification,
      requested_flow_type = params$flow,
      requested_mean_hash = flow_vector_hash(params$flow_mean),
      requested_covariance_hash = flow_matrix_hash(params$flow_cov),
      requested_covariance_dimensions = dim(as.matrix(params$flow_cov)),
      preflight_flow_type = resolved_preflight$flow_type,
      preflight_inner_flow_type = resolved_preflight$inner_flow_type,
      preflight_mean_hash = flow_vector_hash(resolved_preflight$mean),
      preflight_covariance_hash = flow_matrix_hash(resolved_preflight$covariance),
      preflight_agreement = same_mean && same_cov,
      requested_mean = as.numeric(params$flow_mean),
      requested_covariance = unname(as.matrix(params$flow_cov)))
  }

  # Pass arguments to Julia
  for (nm in names(params))
    JuliaCall::julia_assign(nm, params[[nm]])
  JuliaCall::julia_assign("_unc_names", unc_names)
  JuliaCall::julia_assign("free0", free0)
  JuliaCall::julia_assign("stored_velocity0", stored_velocity0)
  JuliaCall::julia_assign("subsampling_anchor_updates", subsampling_anchor_updates)
  JuliaCall::julia_assign("subsampling_warmup", subsampling_warmup)
  JuliaCall::julia_assign("subsampling_use_anchor_bank", subsampling_use_anchor_bank)
  JuliaCall::julia_assign("subsampling_anchor_bank_capacity", subsampling_anchor_bank_capacity)
  JuliaCall::julia_assign("subsampling_anchor_selection_policy",
                          subsampling_anchor_selection_policy)
  JuliaCall::julia_assign("subsampling_position_anchor_bank_spec",
                          subsampling_position_anchor_bank_spec)
  JuliaCall::julia_assign("subsampling_main_anchor_refresh_distance",
                          subsampling_main_anchor_refresh_distance)
  JuliaCall::julia_assign("subsampling_node_source_check_distance",
                          subsampling_node_source_check_distance)
  JuliaCall::julia_assign("subsampling_require_no_main_full_data",
                          subsampling_require_no_main_full_data)
  JuliaCall::julia_assign("subsampling_diagnostic_window",
                          subsampling_diagnostic_window)
  JuliaCall::julia_assign("warmup_adaptation_interval",
                          warmup_adaptation_interval)
  JuliaCall::julia_assign("boomerang_variance_floor_fraction",
                          boomerang_variance_floor_fraction)
  JuliaCall::julia_assign("refresh_rate", as.numeric(refresh_rate))
  JuliaCall::julia_assign("support_boundary_mode", support_boundary$mode)
  JuliaCall::julia_assign("support_boundary_max_bisection_steps", support_boundary$max_bisection_steps)
  JuliaCall::julia_assign("support_boundary_time_rtol", support_boundary$time_rtol)
  JuliaCall::julia_assign("support_boundary_time_atol", support_boundary$time_atol)
  JuliaCall::julia_assign("support_boundary_clip_fraction", support_boundary$clip_fraction)
  JuliaCall::julia_assign("support_boundary_max_refresh_attempts", support_boundary$max_refresh_attempts)
  JuliaCall::julia_assign("support_boundary_refresh_probe_time", support_boundary$refresh_probe_time)
  JuliaCall::julia_assign("support_boundary_min_safe_time", support_boundary$min_safe_time)
  JuliaCall::julia_assign("streaming_trace_directory",
    if (identical(trace_storage, "streaming"))
      normalizePath(streaming_trace_directory, mustWork = FALSE) else NULL)
  JuliaCall::julia_assign("streaming_trace_buffer_events",
                          streaming_trace_buffer_events)

  bridge_profile_mark("argument_transfer_and_validation_seconds")
  if (is.null(subsampling)) {
    result <- .pdmpsamplers_julia_eval("PDMPSamplersRBridge.r_pdmp_stan(
      _pdmp_model, x0, flow, algorithm, flow_mean, flow_cov;
      theta0 = theta0,
      free0 = free0, stored_velocity0 = stored_velocity0,
      c0 = c0, grid_n = grid_n, grid_t_max = grid_t_max,
      use_fd_hvp = _curvature_backend == \"finite_difference\",
      curvature_backend = _curvature_backend,
      post_warmup_simplify = post_warmup_simplify,
      grid_bound = grid_bound,
      grid_curvature_bound = grid_curvature_bound,
      linear_area_threshold = linear_area_threshold,
      linear_min_area_gain = linear_min_area_gain,
      lazy_low_tightness_threshold = lazy_low_tightness_threshold,
      lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
      lazy_max_rejections = lazy_max_rejections,
      t0 = t0, T = T, t_warmup = t_warmup,
      warmup_adaptation_interval = warmup_adaptation_interval,
      boomerang_variance_floor_fraction = boomerang_variance_floor_fraction,
      sticky = sticky, can_stick = can_stick,
      model_prior = model_prior, parameter_prior = parameter_prior,
      slab_prior = slab_prior,
      unc_names = _unc_names,
      show_progress = show_progress, n_chains = n_chains, threaded = threaded,
      seed = seed,
      adaptive_scheme = adaptive_scheme,
      refresh_rate = refresh_rate,
      support_boundary_mode = support_boundary_mode,
      support_boundary_max_bisection_steps = support_boundary_max_bisection_steps,
      support_boundary_time_rtol = support_boundary_time_rtol,
      support_boundary_time_atol = support_boundary_time_atol,
      support_boundary_clip_fraction = support_boundary_clip_fraction,
      support_boundary_max_refresh_attempts = support_boundary_max_refresh_attempts,
      support_boundary_refresh_probe_time = support_boundary_refresh_probe_time,
      support_boundary_min_safe_time = support_boundary_min_safe_time,
      streaming_trace_directory = streaming_trace_directory,
      streaming_trace_buffer_events = streaming_trace_buffer_events
    );")
  } else {
    result <- .pdmpsamplers_julia_eval("PDMPSamplersRBridge.r_pdmp_stan_subsampling(
      _prepared_stan_subsampling, _subsampling, x0, flow, algorithm,
      flow_mean, flow_cov;
      theta0 = theta0,
      free0 = free0, stored_velocity0 = stored_velocity0,
      c0 = c0, grid_n = grid_n, grid_t_max = grid_t_max,
      use_fd_hvp = _curvature_backend == \"finite_difference\",
      curvature_backend = _curvature_backend,
      post_warmup_simplify = post_warmup_simplify,
      grid_bound = grid_bound,
      grid_curvature_bound = grid_curvature_bound,
      linear_area_threshold = linear_area_threshold,
      linear_min_area_gain = linear_min_area_gain,
      lazy_low_tightness_threshold = lazy_low_tightness_threshold,
      lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
      lazy_max_rejections = lazy_max_rejections,
      t0 = t0, T = T, t_warmup = t_warmup,
      warmup_adaptation_interval = warmup_adaptation_interval,
      boomerang_variance_floor_fraction = boomerang_variance_floor_fraction,
      sticky = sticky, can_stick = can_stick,
      model_prior = model_prior, parameter_prior = parameter_prior,
      slab_prior = slab_prior,
      show_progress = show_progress, n_chains = n_chains,
      threaded = threaded, seed = seed,
      adaptive_scheme = adaptive_scheme,
      refresh_rate = refresh_rate,
      warmup_mode = subsampling_warmup,
      n_anchor_updates = subsampling_anchor_updates,
      use_anchor_bank = subsampling_use_anchor_bank,
      anchor_bank_capacity = subsampling_anchor_bank_capacity,
      anchor_selection_policy = subsampling_anchor_selection_policy,
      position_anchor_bank_spec = subsampling_position_anchor_bank_spec,
      main_anchor_refresh_distance = subsampling_main_anchor_refresh_distance,
      node_source_check_distance = subsampling_node_source_check_distance,
      require_no_main_full_data = subsampling_require_no_main_full_data,
      diagnostic_window_width = subsampling_diagnostic_window,
      support_boundary_mode = support_boundary_mode,
      support_boundary_max_bisection_steps = support_boundary_max_bisection_steps,
      support_boundary_time_rtol = support_boundary_time_rtol,
      support_boundary_time_atol = support_boundary_time_atol,
      support_boundary_clip_fraction = support_boundary_clip_fraction,
      support_boundary_max_refresh_attempts = support_boundary_max_refresh_attempts,
      support_boundary_refresh_probe_time = support_boundary_refresh_probe_time,
      support_boundary_min_safe_time = support_boundary_min_safe_time,
      streaming_trace_directory = streaming_trace_directory,
      streaming_trace_buffer_events = streaming_trace_buffer_events
    );")
  }
  bridge_profile_mark("sampling_pack_and_transfer_seconds")
  bridge_internal_profile <- if (bridge_profile_enabled)
    result$bridge_profile else NULL
  if (is.environment(bridge_internal_profile))
    bridge_internal_profile <- as.list(bridge_internal_profile)
  backend_provenance <- if (!is.null(full_gradient)) {
    gradient_type <- tryCatch(as.character(JuliaCall::julia_eval(
      "string(typeof(_pdmp_model.grad.f))")),
      error = function(e) NA_character_)
    hvp_type <- tryCatch(as.character(JuliaCall::julia_eval(
      "string(typeof(_pdmp_model.hvp))")),
      error = function(e) NA_character_)
    resolved <- if (grepl("OMRFFullGradientProvider", gradient_type,
                           fixed = TRUE)) "native" else "unresolved"
    list(requested_backend = "native", resolved_backend = resolved,
         gradient_provider_type = gradient_type,
         hvp_provider_type = hvp_type)
  } else {
    list(requested_backend = "bridgestan", resolved_backend = "bridgestan",
         gradient_provider_type = NA_character_,
         hvp_provider_type = NA_character_)
  }
  if (is.environment(result)) result <- as.list(result)
  if (verify_flow_specification) {
    resolved_chains <- lapply(seq_len(params$n_chains), function(chain) {
      value <- .pdmpsamplers_julia_call(
        "r_chain_flow_spec", result$chains, chain = as.integer(chain))
      if (is.environment(value)) value <- as.list(value)
      value
    })
    agreement <- .flow_spec_agreement(resolved_chains,
      requested_flow_provenance$requested_mean,
      requested_flow_provenance$requested_covariance)
    mean_agreement <- agreement$mean
    covariance_agreement <- agreement$covariance
    chain_agreement <- agreement$chain
    if (resolved_flow_specification == "frozen" && !all(chain_agreement)) {
      mismatched <- paste(which(!chain_agreement), collapse = ", ")
      cli::cli_abort(c(
        "At least one sampled Julia chain did not retain the requested frozen flow.",
        "i" = "Mismatched chain indices: {mismatched}.",
        "i" = "Refusing to return a frozen-flow result."
      ))
    }
    requested_flow_provenance$resolved_flow_type <- vapply(
      resolved_chains, `[[`, "", "flow_type")
    requested_flow_provenance$resolved_inner_flow_type <- vapply(
      resolved_chains, `[[`, "", "inner_flow_type")
    requested_flow_provenance$resolved_mean_hash <- vapply(
      resolved_chains, function(x) flow_vector_hash(x$mean), "")
    requested_flow_provenance$resolved_covariance_hash <- vapply(
      resolved_chains, function(x) flow_matrix_hash(x$covariance), "")
    requested_flow_provenance$resolved_mean_agreement <- mean_agreement
    requested_flow_provenance$resolved_covariance_agreement <- covariance_agreement
    requested_flow_provenance$resolved_agreement <- all(chain_agreement)
    requested_flow_provenance$verification_passed <-
      resolved_flow_specification == "initial" || all(chain_agreement)
    requested_flow_provenance$n_verified_chains <- length(resolved_chains)
  }
  subsampling_context_counters <- result$subsampling_context_counters
  window_diagnostics <- result$window_diagnostics
  candidate_records <- result$candidate_records
  if (!is.null(subsampling_context_counters)) {
    subsampling_context_counters <- lapply(subsampling_context_counters, function(counter) {
      if (is.environment(counter)) as.list(counter) else counter
    })
  }
  result <- new_pdmp_result(
    chains   = result$chains,
    stats    = result$stats,
    d        = result$d,
    n_chains = result$n_chains
  )
  attr(result, "backend_provenance") <- backend_provenance
  attr(result, "flow_provenance") <- requested_flow_provenance
  if (!is.null(subsampling)) {
    attr(result, "subsampling") <- TRUE
    attr(result, "subsampling_context_counters") <- subsampling_context_counters
    attr(result, "omrf_window_diagnostics") <- window_diagnostics
    attr(result, "omrf_candidate_records") <- candidate_records
  }
  if (materialize) result <- .materialize(result)
  if (!is.null(subsampling) &&
      subsampling_anchor_selection_policy %in% c("position", "node_position")) {
    dataset_hash <- dataset_identity_hash(standata)
    bank_envelope_spec <- subsampling$residual_envelope
    bank_envelope_spec$categorical_residual_cap <- NULL
    envelope_hash <- flow_hash(bank_envelope_spec)
    node_support_order <- function(envelope, parameter_order) {
      seen <- as.integer(envelope$seen)
      P <- length(seen)
      thresholds <- grep(paste0("^", envelope$thresholds, "\\."),
        parameter_order)
      interactions <- grep(paste0("^", envelope$interactions, "\\."),
        parameter_order)
      starts <- cumsum(c(1L, head(seen - 1L, -1L)))
      incident <- vector("list", P)
      edge <- 0L
      if (P > 1L) for (left in seq_len(P - 1L))
        for (right in (left + 1L):P) {
          edge <- edge + 1L
          incident[[left]] <- c(incident[[left]], interactions[[edge]])
          incident[[right]] <- c(incident[[right]], interactions[[edge]])
        }
      lapply(seq_len(P), function(node) c(
        thresholds[seq.int(starts[[node]], length.out = seen[[node]] - 1L)],
        incident[[node]]))
    }
    supplied_bank_specs <- if (is.null(subsampling_position_anchor_bank_spec))
      NULL else if (!is.null(subsampling_position_anchor_bank_spec$version))
        rep(list(subsampling_position_anchor_bank_spec), params$n_chains) else
        subsampling_position_anchor_bank_spec
    bank_specs <- lapply(seq_len(params$n_chains), function(chain) {
      counters <- if (length(subsampling_context_counters) >= chain &&
          is.list(subsampling_context_counters[[chain]]))
        subsampling_context_counters[[chain]] else list()
      endpoint <- if (!is.null(result$terminal_states))
        result$terminal_states[[chain]] else NULL
      supplied <- if (is.null(supplied_bank_specs)) NULL else
        supplied_bank_specs[[chain]]
      .position_bank_result_spec(counters, supplied, endpoint,
        subsampling_anchor_selection_policy, unc_names, dataset_hash,
        envelope_hash, subsampling$residual_envelope, node_support_order)
    })
    attr(result, "position_anchor_bank_spec") <- if (length(bank_specs) == 1L)
      bank_specs[[1L]] else bank_specs
  }
  bridge_profile_mark("output_metadata_seconds")
  if (bridge_profile_enabled) {
    bridge_profile$julia_internal <- bridge_internal_profile
    bridge_profile$total_seconds <- proc.time()[["elapsed"]] -
      bridge_profile_started
    attr(result, "bridge_phase_profile") <- bridge_profile
  }
  result
}
