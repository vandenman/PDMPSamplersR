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

# Internal function to validate and prepare common PDMP parameters
validate_pdmp_params <- function(d, flow, algorithm, T, t0 = 0.0, t_warmup = 0.0,
                                flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                                x0 = NULL, theta0 = NULL, show_progress = TRUE,
                                sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                                slab_prior = NULL,
                                grid_n = 30, grid_t_max = 2.0,
                                post_warmup_simplify = FALSE,
                                grid_bound = "constant",
                                linear_area_threshold = 0.95,
                                linear_min_area_gain = 0.0,
                                n_chains = 1L, threaded = FALSE, seed = NULL,
                                adaptive_scheme = "diagonal",
                                unc_names = NULL,
                                lazy_low_tightness_threshold = 0.1,
                                lazy_max_low_tightness_rejections = 3L,
                                lazy_max_rejections = 0L) {

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
  grid_like_algorithms <- c("GridThinningStrategy", "PositiveVariationGridThinningStrategy",
                            "VectorVariationThinningStrategy")
  algorithm <- match.arg(algorithm, c("ThinningStrategy", grid_like_algorithms, "RootsPoissonStrategy"))
  adaptive_scheme <- match.arg(adaptive_scheme, c("diagonal", "fullrank"))
  grid_bound <- match.arg(grid_bound, c("constant", "flat", "linear", "auto"))

  # AdaptiveBoomerang requires GridThinningStrategy and a warmup period
  if (flow == "AdaptiveBoomerang") {
    if (!algorithm %in% grid_like_algorithms)
      cli::cli_abort("{.val AdaptiveBoomerang} requires a grid-like algorithm.")
    if (t_warmup == 0) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for {.val AdaptiveBoomerang}.")
    }
  }

  # PreconditionedZigZag and PreconditionedBPS require GridThinningStrategy and a warmup period
  if (flow %in% c("PreconditionedZigZag", "PreconditionedBPS")) {
    if (!algorithm %in% grid_like_algorithms)
      cli::cli_abort("{.val {flow}} requires a grid-like algorithm.")
    if (t_warmup == 0) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for {.val {flow}}.")
    }
  }

  validate_type(c0, type = "double", n = 1, positive = TRUE)
  validate_type(show_progress, type = "logical", n = 1)

  n_chains <- cast_integer(n_chains, n = 1)
  validate_type(n_chains, type = "integer", n = 1, positive = TRUE)
  validate_type(threaded, type = "logical", n = 1)

  if (threaded && rlang::is_false(.pdmpsamplers_julia_eval("r_threading_available()"))) {
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
      if (!flow %in% c("ZigZag", "BouncyParticle")) {
        cli::cli_abort("Dependent {.arg slab_prior} sticky sampling currently supports only ZigZag and BouncyParticle flows.")
      }
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
    idx <- match(coef, unc_names)
    if (anyNA(idx)) {
      missing <- coef[is.na(idx)]
      cli::cli_abort(c(
        "{.arg coef} contains names not found among unconstrained parameter names.",
        "x" = "Unknown names: {.val {missing}}."
      ))
    }
    return(idx)
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
    idx <- match(value, unc_names)
    if (anyNA(idx)) {
      missing <- value[is.na(idx)]
      cli::cli_abort(c(
        "{.arg {arg}} contains names not found among unconstrained parameter names.",
        "x" = "Unknown names: {.val {missing}}."
      ))
    }
    return(idx)
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
  } else if (slab_prior$type == "global_logscale_exchangeable_gaussian") {
    logscale_idx <- .slab_state_indices(slab_prior$logscale, d, unc_names, arg = "logscale")
    if (length(logscale_idx) != 1L) {
      cli::cli_abort("Argument {.arg logscale} must have length 1.")
    }
    beta_idx <- if (is.null(coef_idx)) which(can_stick) else coef_idx
    if (logscale_idx %in% beta_idx) {
      cli::cli_abort("Slab {.arg logscale} coordinate must be disjoint from slab {.arg coef} coordinates.")
    }
    if (can_stick[logscale_idx]) {
      cli::cli_abort("Slab {.arg logscale} coordinate must be non-stickable.")
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
                                           model_prior, parameter_prior,
                                           subsample = NULL) {
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
  if (!flow %in% c("ZigZag", "BouncyParticle")) {
    cli::cli_abort("Dependent {.arg slab_prior} sticky sampling currently supports only ZigZag and BouncyParticle flows.")
  }
  if (algorithm != "GridThinningStrategy") {
    cli::cli_abort("Dependent {.arg slab_prior} sticky sampling currently requires {.val GridThinningStrategy}.")
  }
  if (!is.null(subsample)) {
    cli::cli_abort("Dependent {.arg slab_prior} is not yet supported together with {.arg subsample}.")
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

.standata_n_obs <- function(standata, standata_path, subsample) {
  if (!is.null(subsample$N)) {
    n <- cast_integer(subsample$N, n = 1)
    validate_type(n, type = "integer", n = 1, positive = TRUE)
    return(n)
  }

  if (is.list(standata) && !is.null(standata$N)) {
    n <- cast_integer(standata$N, n = 1)
    validate_type(n, type = "integer", n = 1, positive = TRUE)
    return(n)
  }

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    data <- tryCatch(
      jsonlite::read_json(standata_path, simplifyVector = TRUE),
      error = function(e) NULL
    )
    if (is.list(data) && !is.null(data$N)) {
      n <- cast_integer(data$N, n = 1)
      validate_type(n, type = "integer", n = 1, positive = TRUE)
      return(n)
    }
  }

  cli::cli_abort(c(
    "Could not determine the number of observations for subsampling.",
    "i" = "Provide {.field N} in {.arg standata} or set {.code subsample$N}."
  ))
}

.write_or_validate_standata_path <- function(data, arg = "subsample$prior_standata") {
  if (is.list(data)) {
    path <- tempfile(fileext = ".json")
    write_stan_json(data, path)
    return(list(path = path, temporary = TRUE))
  }

  validate_type(data, type = "character", n = 1)
  if (!file.exists(data))
    cli::cli_abort("{.arg {arg}} file not found: {.path {data}}")
  if (!grepl("\\.json$", data))
    cli::cli_abort("{.arg {arg}} should be a JSON file path or a named list.")
  list(path = data, temporary = FALSE)
}

validate_stan_subsample <- function(subsample, standata, standata_path, path_to_stanmodel) {
  if (is.null(subsample)) return(NULL)
  if (!is.list(subsample))
    cli::cli_abort("{.arg subsample} must be NULL or a named list.")

  if (is.null(subsample$size))
    cli::cli_abort("{.arg subsample$size} must be provided.")
  N <- .standata_n_obs(standata, standata_path, subsample)
  size <- cast_integer(subsample$size, n = 1)
  validate_type(size, type = "integer", n = 1, positive = TRUE)
  if (size >= N)
    cli::cli_abort("{.arg subsample$size} ({size}) must be between 1 and {.code N - 1} ({N - 1L}).")

  prior <- subsample$prior_standata
  if (is.null(prior)) prior <- subsample$prior_data
  if (is.null(prior)) prior <- subsample$prior_standata_path
  if (is.null(prior))
    cli::cli_abort("{.arg subsample$prior_standata} must be provided as a named list or JSON path.")
  prior_info <- .write_or_validate_standata_path(prior, "subsample$prior_standata")

  stanmodel_sub <- subsample$path_to_stanmodel
  if (is.null(stanmodel_sub)) stanmodel_sub <- subsample$stan_file_ext
  if (is.null(stanmodel_sub)) stanmodel_sub <- subsample$stanmodel
  if (is.null(stanmodel_sub)) stanmodel_sub <- path_to_stanmodel
  validate_type(stanmodel_sub, type = "character", n = 1)
  if (!file.exists(stanmodel_sub))
    cli::cli_abort("Subsampled Stan model file not found: {.path {stanmodel_sub}}")
  if (!grepl("\\.(so|dll|dylib|stan)$", stanmodel_sub))
    cli::cli_abort("{.arg subsample$path_to_stanmodel} should point to a Stan model or compiled Stan library.")

  hpp <- subsample$hpp_path
  if (is.null(hpp)) hpp <- subsample$external_hpp
  if (is.null(hpp)) hpp <- hpp_path()
  validate_type(hpp, type = "character", n = 1)
  if (!file.exists(hpp))
    cli::cli_abort("Subsampled Stan external header not found: {.path {hpp}}")

  n_anchor_updates <- if (is.null(subsample$n_anchor_updates)) 10L else cast_integer(subsample$n_anchor_updates, n = 1)
  validate_type(n_anchor_updates, type = "integer", n = 1)
  if (n_anchor_updates < 0)
    cli::cli_abort("{.arg subsample$n_anchor_updates} must be non-negative.")

  hvp_mode <- if (is.null(subsample$hvp_mode)) "scaled" else subsample$hvp_mode
  validate_type(hvp_mode, type = "character", n = 1)
  if (!hvp_mode %in% c("scaled", "none"))
    cli::cli_abort("{.arg subsample$hvp_mode} must be one of {.val scaled} or {.val none}.")

  bool_option <- function(name, default) {
    value <- subsample[[name]]
    if (is.null(value)) value <- default
    validate_type(value, type = "logical", n = 1)
    value
  }
  double_option <- function(name, default, non_negative = TRUE) {
    value <- subsample[[name]]
    if (is.null(value)) value <- default
    validate_type(value, type = "double", n = 1)
    if (non_negative && value < 0)
      cli::cli_abort("{.arg {paste0('subsample$', name)}} must be non-negative.")
    value
  }
  integer_option <- function(name, default, positive = TRUE) {
    value <- subsample[[name]]
    if (is.null(value)) value <- default
    value <- cast_integer(value, n = 1)
    validate_type(value, type = "integer", n = 1)
    if (positive && value <= 0)
      cli::cli_abort("{.arg {paste0('subsample$', name)}} must be positive.")
    value
  }

  list(
    N = N,
    size = size,
    prior_path = prior_info$path,
    prior_temporary = prior_info$temporary,
    path_to_stanmodel = stanmodel_sub,
    hpp_path = hpp,
    n_anchor_updates = n_anchor_updates,
    hvp_mode = hvp_mode,
    use_hcv = bool_option("use_hcv", FALSE),
    use_anchor_bank = bool_option("use_anchor_bank", FALSE),
    use_fd_hvp = bool_option("use_fd_hvp", FALSE),
    compute_lp = bool_option("compute_lp", FALSE),
    resample_dt = double_option("resample_dt", 0.0),
    discretize_dt = double_option("discretize_dt", 0.0),
    bank_capacity = integer_option("bank_capacity", 20L),
    use_fd_hcv = bool_option("use_fd_hcv", FALSE),
    output_csv = if (is.null(subsample$output_csv)) tempfile(fileext = ".csv") else {
      validate_type(subsample$output_csv, type = "character", n = 1)
      subsample$output_csv
    }
  )
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
#'   during warmup and requires a grid-like algorithm.
#'   The `"PreconditionedZigZag"` and `"PreconditionedBPS"` flows learn a
#'   diagonal preconditioner during warmup and also require
#'   a grid-like algorithm.
#' @param algorithm Character string specifying the algorithm. One of
#'   "ThinningStrategy", "GridThinningStrategy",
#'   "PositiveVariationGridThinningStrategy", "VectorVariationThinningStrategy",
#'   or "RootsPoissonStrategy".
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
#'   \code{parameter_prior}. The custom-gradient path requires
#'   \code{prior_grad} so the dependent target can subtract the slab component
#'   before adding the active slab gradient.
#' @param prior_grad Optional function computing the NEGATIVE gradient of the
#'   slab prior component to subtract in the dependent-slab target. Do not
#'   include nuisance-prior terms here; nuisance coordinates should remain in
#'   the posterior gradient \code{f}. Required when \code{slab_prior} is
#'   supplied to \code{pdmp_sample()}.
#' @param grid_n Integer, number of grid points for GridThinningStrategy (default: 30).
#' @param grid_t_max Numeric, maximum time for grid in GridThinningStrategy (default: 2.0).
#' @param post_warmup_simplify Logical. If \code{TRUE}, the grid-thinning
#'   strategy may switch to a constant bound after warmup, reducing
#'   gradient calls during the main sampling phase.
#' @param grid_bound Character string selecting the GridThinning proposal
#'   envelope. \code{"constant"} preserves the historical behavior.
#'   \code{"flat"}, \code{"linear"}, and \code{"auto"} expose signed-rate
#'   affine envelope machinery when supported by the chosen flow and derivative
#'   provider.
#' @param linear_area_threshold Numeric gate for \code{grid_bound = "linear"}
#'   or \code{"auto"}; affine cells are skipped when their relative area gain is
#'   too small.
#' @param linear_min_area_gain Numeric absolute area-gain gate for affine cells.
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
                                      "PositiveVariationGridThinningStrategy", "VectorVariationThinningStrategy",
                                      "RootsPoissonStrategy"),
                        T = 50000, t0 = 0.0, t_warmup = 0.0,
                        flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                        x0 = NULL, theta0 = NULL,
                        hessian = NULL, prior_grad = NULL,
                        sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                        slab_prior = NULL,
                        grid_n = 30, grid_t_max = 2.0,
                        use_fd_hvp = FALSE,
                        post_warmup_simplify = FALSE,
                        grid_bound = c("constant", "flat", "linear", "auto"),
                        linear_area_threshold = 0.95,
                        linear_min_area_gain = 0.0,
                        lazy_low_tightness_threshold = 0.1,
                        lazy_max_low_tightness_rejections = 3L,
                        lazy_max_rejections = 0L,
                        show_progress = TRUE,
                        n_chains = 1L, threaded = FALSE, seed = NULL,
                        adaptive_scheme = c("diagonal", "fullrank"),
                        materialize = TRUE,
                        support_boundary = support_boundary_control()) {

  # Validate function argument (fail fast before Julia setup)
  if (!rlang::is_function(f)) {
    cli::cli_abort(c(
      "{.arg f} must be a function.",
      "x" = "You've supplied an object of {.cls {class(f)}}."
    ))
  }

  # Use common validation function
  params <- validate_pdmp_params(d, flow, algorithm, T, t0, t_warmup, flow_mean, flow_cov,
                                c0, x0, theta0, show_progress,
                                sticky, can_stick, model_prior, parameter_prior, slab_prior,
                                grid_n, grid_t_max, post_warmup_simplify,
                                grid_bound, linear_area_threshold, linear_min_area_gain,
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
  if (!is.null(slab_prior)) {
    if (!rlang::is_function(prior_grad)) {
      cli::cli_abort("Argument {.arg prior_grad} must be provided as a slab-only negative-gradient function when {.arg slab_prior} is supplied.")
    }
    tryCatch({
      test_output <- prior_grad(params$x0)
      if (!is.numeric(test_output) || length(test_output) != params$d) {
        cli::cli_abort("{.arg prior_grad} must return a numeric vector of length {params$d}, but returned length {length(test_output)}.")
      }
    }, error = function(e) {
      cli::cli_abort(c(
        "{.arg prior_grad} failed on test input.",
        "x" = conditionMessage(e)
      ))
    })
  }

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
  JuliaCall::julia_assign("prior_grad_r", prior_grad)
  JuliaCall::julia_command("prior_grad!(out, x) = out .= prior_grad_r(x);")
  JuliaCall::julia_assign("hessian_f", hessian)
  JuliaCall::julia_assign("support_boundary_mode", support_boundary$mode)
  JuliaCall::julia_assign("support_boundary_max_bisection_steps", support_boundary$max_bisection_steps)
  JuliaCall::julia_assign("support_boundary_time_rtol", support_boundary$time_rtol)
  JuliaCall::julia_assign("support_boundary_time_atol", support_boundary$time_atol)
  JuliaCall::julia_assign("support_boundary_clip_fraction", support_boundary$clip_fraction)
  JuliaCall::julia_assign("support_boundary_max_refresh_attempts", support_boundary$max_refresh_attempts)
  JuliaCall::julia_assign("support_boundary_refresh_probe_time", support_boundary$refresh_probe_time)
  JuliaCall::julia_assign("support_boundary_min_safe_time", support_boundary$min_safe_time)

  result <- .pdmpsamplers_julia_eval("r_pdmp_custom(
    grad!, d, x0, flow, algorithm, flow_mean, flow_cov;
    c0 = c0, grid_n = grid_n, grid_t_max = grid_t_max,
    post_warmup_simplify = post_warmup_simplify,
    grid_bound = grid_bound,
    linear_area_threshold = linear_area_threshold,
    linear_min_area_gain = linear_min_area_gain,
    t0 = t0, T = T, t_warmup = t_warmup,
    hessian = hessian_f,
    sticky = sticky, can_stick = can_stick,
    model_prior = model_prior, parameter_prior = parameter_prior,
    slab_prior = slab_prior,
    prior_grad! = prior_grad!,
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
    support_boundary_min_safe_time = support_boundary_min_safe_time
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
#' @param prior_stanmodel Optional Stan model path for the base prior target
#'   used by dependent slabs. Defaults to \code{path_to_stanmodel} when
#'   \code{prior_standata} is supplied.
#' @param prior_standata Optional prior-only Stan data path or named list for
#'   the base prior target. Required when \code{slab_prior} is supplied.
#' @param subsample NULL (default) for full-data sampling, or a named list with
#'   at least `size` and `prior_standata`. Optional entries include
#'   `path_to_stanmodel` (the external-C++ subsampled Stan model; defaults to
#'   `path_to_stanmodel`), `hpp_path`, `n_anchor_updates`, `hvp_mode`,
#'   `use_hcv`, `use_anchor_bank`, `use_fd_hvp`, `compute_lp`, and
#'   `resample_dt`.
#' @inheritParams pdmp_sample
#'
#' @return A \code{pdmp_result} object. Use \code{mean}, \code{var},
#'   \code{quantile}, etc. for continuous-time estimators, or
#'   \code{discretize} to obtain a sample matrix.
#'
#' @export
pdmp_sample_from_stanmodel <- function(path_to_stanmodel, standata,
                        prior_stanmodel = NULL, prior_standata = NULL,
                        flow = c("ZigZag", "BouncyParticle", "Boomerang", "AdaptiveBoomerang", "PreconditionedZigZag", "PreconditionedBPS"),
                        algorithm = c("ThinningStrategy", "GridThinningStrategy",
                                      "PositiveVariationGridThinningStrategy", "VectorVariationThinningStrategy",
                                      "RootsPoissonStrategy"),
                        T = 50000, t0 = 0.0, t_warmup = 0.0,
                        flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                        x0 = NULL, theta0 = NULL,
                        sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                        slab_prior = NULL,
                        grid_n = 30, grid_t_max = 2.0,
                        use_fd_hvp = FALSE,
                        post_warmup_simplify = FALSE,
                        grid_bound = c("constant", "flat", "linear", "auto"),
                        linear_area_threshold = 0.95,
                        linear_min_area_gain = 0.0,
                        lazy_low_tightness_threshold = 0.1,
                        lazy_max_low_tightness_rejections = 3L,
                        lazy_max_rejections = 0L,
                        show_progress = TRUE,
                        n_chains = 1L, threaded = FALSE, seed = NULL,
                        adaptive_scheme = c("diagonal", "fullrank"),
                        materialize = TRUE,
                        support_boundary = support_boundary_control(),
                        subsample = NULL) {

  # Validate file paths on R side before setting up Julia
  validate_type(path_to_stanmodel, type = "character", n = 1)
  flow <- match.arg(flow)
  algorithm <- match.arg(algorithm)
  adaptive_scheme <- match.arg(adaptive_scheme)

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
  .validate_dependent_slab_early(
    slab_prior, sticky, flow, algorithm, model_prior, parameter_prior,
    subsample = subsample
  )
  if (!is.null(slab_prior)) {
    cli::cli_abort("Dependent {.arg slab_prior} is temporarily gated for Stan models until target composition subtracts only the slab component or adds back nuisance priors.")
  }
  if (!is.null(slab_prior) && is.null(prior_standata_path)) {
    cli::cli_abort("Argument {.arg prior_standata} is required when {.arg slab_prior} is supplied.")
  }
  if (is.null(prior_stanmodel)) {
    prior_stanmodel <- path_to_stanmodel
  } else {
    validate_type(prior_stanmodel, type = "character", n = 1)
  }

  if (!file.exists(path_to_stanmodel))
    cli::cli_abort("Stan model file not found: {.path {path_to_stanmodel}}")
  if (!file.exists(standata_path))
    cli::cli_abort("Stan data file not found: {.path {standata_path}}")
  if (!is.null(slab_prior)) {
    if (!file.exists(prior_stanmodel))
      cli::cli_abort("Prior Stan model file not found: {.path {prior_stanmodel}}")
    if (!file.exists(prior_standata_path))
      cli::cli_abort("Prior Stan data file not found: {.path {prior_standata_path}}")
  }
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
  if (!is.null(slab_prior)) {
    if (!grepl("\\.(so|dll|dylib|stan)$", prior_stanmodel)) {
      cli::cli_abort("{.arg prior_stanmodel} should point to a Stan model or compiled Stan library.")
    }
    if (!grepl("\\.json$", prior_standata_path)) {
      cli::cli_abort("{.arg prior_standata} should be a JSON file path or a list that can be written to JSON.")
    }
  }

  subsample <- validate_stan_subsample(subsample, standata, standata_path, path_to_stanmodel)
  if (!is.null(subsample) && isTRUE(subsample$prior_temporary))
    on.exit(unlink(subsample$prior_path), add = TRUE)

  check_for_julia_setup()

  # Normalize paths to absolute
  path_to_stanmodel <- normalizePath(path_to_stanmodel, mustWork = TRUE)
  standata_path     <- normalizePath(standata_path,     mustWork = TRUE)
  if (!is.null(slab_prior)) {
    prior_stanmodel <- normalizePath(prior_stanmodel, mustWork = TRUE)
    prior_standata_path <- normalizePath(prior_standata_path, mustWork = TRUE)
  }

  support_boundary <- validate_support_boundary_control(support_boundary)
  validate_type(use_fd_hvp, type = "logical", n = 1)

  JuliaCall::julia_assign("_path_to_stan_model", path_to_stanmodel)
  JuliaCall::julia_assign("_path_to_stan_data",  standata_path)
  JuliaCall::julia_assign("_path_to_prior_stan_model", prior_stanmodel)
  JuliaCall::julia_assign("_path_to_prior_stan_data_full", prior_standata_path)

  if (is.null(subsample)) {
    # Create PDMPModel in Julia and get dimension for the ordinary full-data path.
    JuliaCall::julia_assign("_use_fd_hvp", isTRUE(use_fd_hvp))
    JuliaCall::julia_command("_stan_model = BridgeStan.StanModel(_path_to_stan_model, _path_to_stan_data; warn=false); _pdmp_model = PDMPModel(_stan_model; hvp = !_use_fd_hvp);")
    d <- JuliaCall::julia_eval("_pdmp_model.d")
    unc_names <- character(0)
    if (!is.null(slab_prior)) {
      JuliaCall::julia_command("_prior_stan_model = BridgeStan.StanModel(_path_to_prior_stan_model, _path_to_prior_stan_data_full; warn=false); _prior_pdmp_model = PDMPModel(_prior_stan_model);")
      d_prior <- JuliaCall::julia_eval("_prior_pdmp_model.d")
      if (!identical(as.integer(d), as.integer(d_prior))) {
        cli::cli_abort("Posterior and prior Stan models must have the same unconstrained dimension.")
      }
      unc_names <- JuliaCall::julia_eval("BridgeStan.param_unc_names(_stan_model)")
      prior_unc_names <- JuliaCall::julia_eval("BridgeStan.param_unc_names(_prior_stan_model)")
      if (!identical(unc_names, prior_unc_names)) {
        cli::cli_abort("Posterior and prior Stan models must have identical unconstrained parameter names.")
      }
    } else {
      JuliaCall::julia_command("_prior_pdmp_model = nothing;")
    }
  } else {
    subsample$path_to_stanmodel <- normalizePath(subsample$path_to_stanmodel, mustWork = TRUE)
    subsample$prior_path <- normalizePath(subsample$prior_path, mustWork = TRUE)
    subsample$hpp_path <- normalizePath(subsample$hpp_path, mustWork = TRUE)
    d <- .pdmpsamplers_julia_call(
      "r_stan_param_unc_num_with_header",
      path_to_stanmodel,
      standata_path,
      subsample$hpp_path
    )
    unc_names <- character(0)
  }

  # Use common validation function
  params <- validate_pdmp_params(d, flow, algorithm, T, t0, t_warmup, flow_mean, flow_cov,
                                 c0, x0, theta0, show_progress,
                                 sticky, can_stick, model_prior, parameter_prior, slab_prior,
                                 grid_n, grid_t_max, post_warmup_simplify,
                                 grid_bound, linear_area_threshold, linear_min_area_gain,
                                 n_chains, threaded, seed,
                                 adaptive_scheme = adaptive_scheme,
                                 unc_names = unc_names,
                                 lazy_low_tightness_threshold = lazy_low_tightness_threshold,
                                 lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
                                 lazy_max_rejections = lazy_max_rejections)

  if (isTRUE(params$threaded) && params$n_chains > 1L) {
    cli::cli_warn(c(
      "Stan-backed PDMP sampling uses BridgeStan gradients, which are serialized across Julia threads to avoid Stan Math autodiff memory corruption.",
      "i" = "This prevents segmentation faults but may limit parallel-chain speedups.",
      "i" = "For true parallel speedups with Stan-backed models, use separate R/Julia processes or a Julia-native thread-safe gradient implementation."
    ))
  }

  # Pass arguments to Julia
  for (nm in names(params))
    JuliaCall::julia_assign(nm, params[[nm]])
  JuliaCall::julia_assign("_unc_names", unc_names)
  JuliaCall::julia_assign("support_boundary_mode", support_boundary$mode)
  JuliaCall::julia_assign("support_boundary_max_bisection_steps", support_boundary$max_bisection_steps)
  JuliaCall::julia_assign("support_boundary_time_rtol", support_boundary$time_rtol)
  JuliaCall::julia_assign("support_boundary_time_atol", support_boundary$time_atol)
  JuliaCall::julia_assign("support_boundary_clip_fraction", support_boundary$clip_fraction)
  JuliaCall::julia_assign("support_boundary_max_refresh_attempts", support_boundary$max_refresh_attempts)
  JuliaCall::julia_assign("support_boundary_refresh_probe_time", support_boundary$refresh_probe_time)
  JuliaCall::julia_assign("support_boundary_min_safe_time", support_boundary$min_safe_time)

  if (is.null(subsample)) {
    result <- .pdmpsamplers_julia_eval("r_pdmp_stan(
      _pdmp_model, x0, flow, algorithm, flow_mean, flow_cov;
      c0 = c0, grid_n = grid_n, grid_t_max = grid_t_max,
      post_warmup_simplify = post_warmup_simplify,
      grid_bound = grid_bound,
      linear_area_threshold = linear_area_threshold,
      linear_min_area_gain = linear_min_area_gain,
      lazy_low_tightness_threshold = lazy_low_tightness_threshold,
      lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
      lazy_max_rejections = lazy_max_rejections,
      t0 = t0, T = T, t_warmup = t_warmup,
      sticky = sticky, can_stick = can_stick,
      model_prior = model_prior, parameter_prior = parameter_prior,
      slab_prior = slab_prior,
      prior_model = _prior_pdmp_model,
      unc_names = _unc_names,
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
      support_boundary_min_safe_time = support_boundary_min_safe_time
    );")
  } else {
    JuliaCall::julia_assign("_path_to_subsampled_stan_model", subsample$path_to_stanmodel)
    JuliaCall::julia_assign("_path_to_prior_stan_data", subsample$prior_path)
    JuliaCall::julia_assign("_subsample_hpp_path", subsample$hpp_path)
    JuliaCall::julia_assign("_subsample_N", as.integer(subsample$N))
    JuliaCall::julia_assign("_subsample_size", as.integer(subsample$size))
    JuliaCall::julia_assign("_subsample_output_csv", subsample$output_csv)
    JuliaCall::julia_assign("_subsample_n_anchor_updates", as.integer(subsample$n_anchor_updates))
    JuliaCall::julia_assign("_subsample_hvp_mode", subsample$hvp_mode)
    JuliaCall::julia_assign("_subsample_use_hcv", subsample$use_hcv)
    JuliaCall::julia_assign("_subsample_use_anchor_bank", subsample$use_anchor_bank)
    JuliaCall::julia_assign("_subsample_use_fd_hvp", subsample$use_fd_hvp)
    JuliaCall::julia_assign("_subsample_compute_lp", subsample$compute_lp)
    JuliaCall::julia_assign("_subsample_resample_dt", subsample$resample_dt)
    JuliaCall::julia_assign("_subsample_discretize_dt", subsample$discretize_dt)
    JuliaCall::julia_assign("_subsample_bank_capacity", as.integer(subsample$bank_capacity))
    JuliaCall::julia_assign("_subsample_use_fd_hcv", subsample$use_fd_hcv)

    result <- .pdmpsamplers_julia_eval("r_pdmp_brms_subsampled(
      _path_to_stan_model, _path_to_subsampled_stan_model, _subsample_hpp_path,
      _path_to_stan_data, _path_to_prior_stan_data,
      _subsample_N, _subsample_size,
      flow, algorithm, flow_mean, flow_cov, _subsample_output_csv;
      c0 = c0, grid_n = grid_n, grid_t_max = grid_t_max,
      grid_bound = grid_bound,
      linear_area_threshold = linear_area_threshold,
      linear_min_area_gain = linear_min_area_gain,
      t0 = t0, T = T, t_warmup = t_warmup,
      n_anchor_updates = _subsample_n_anchor_updates,
      adaptive_scheme = adaptive_scheme,
      discretize_dt = _subsample_discretize_dt,
      show_progress = show_progress,
      n_chains = n_chains, threaded = threaded, seed = seed,
      compute_lp = _subsample_compute_lp,
      resample_dt = _subsample_resample_dt,
      hvp_mode = _subsample_hvp_mode,
      use_hcv = _subsample_use_hcv,
      use_anchor_bank = _subsample_use_anchor_bank,
      bank_capacity = _subsample_bank_capacity,
      use_fd_hvp = _subsample_use_fd_hvp,
      post_warmup_simplify = post_warmup_simplify,
      use_fd_hcv = _subsample_use_fd_hcv,
      sticky = sticky, can_stick = can_stick,
      model_prior = model_prior, parameter_prior = parameter_prior,
      slab_prior = slab_prior,
      unc_names = _unc_names
    );")
  }
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
