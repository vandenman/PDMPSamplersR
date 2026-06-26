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
                                grid_n = 30, grid_t_max = 2.0,
                                post_warmup_simplify = FALSE,
                                n_chains = 1L, threaded = FALSE, seed = NULL,
                                adaptive_scheme = "diagonal") {

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
  algorithm <- match.arg(algorithm, c("ThinningStrategy", "GridThinningStrategy", "RootsPoissonStrategy"))
  adaptive_scheme <- match.arg(adaptive_scheme, c("diagonal", "fullrank"))

  # AdaptiveBoomerang requires GridThinningStrategy and a warmup period
  if (flow == "AdaptiveBoomerang") {
    if (algorithm != "GridThinningStrategy")
      cli::cli_abort("{.val AdaptiveBoomerang} requires {.val GridThinningStrategy} as the algorithm.")
    if (t_warmup == 0) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for {.val AdaptiveBoomerang}.")
    }
  }

  # PreconditionedZigZag and PreconditionedBPS require GridThinningStrategy and a warmup period
  if (flow %in% c("PreconditionedZigZag", "PreconditionedBPS")) {
    if (algorithm != "GridThinningStrategy")
      cli::cli_abort("{.val {flow}} requires {.val GridThinningStrategy} as the algorithm.")
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
  if (sticky) {
    if (is.null(can_stick)) {
      can_stick <- rep(FALSE, d)
    } else {
      validate_type(can_stick, type = "logical", n = d)
    }
    if (is.null(model_prior) || !(is.bernoulli(model_prior) || is.betabernoulli(model_prior))) {
      cli::cli_abort("Argument {.arg model_prior} must be provided when {.arg sticky} is {.code TRUE}. It should be an object of class {.cls bernoulli} or {.cls beta-bernoulli}.")
    } else if (is.bernoulli(model_prior)) {
      if (length(model_prior$prob) == 1) {
        model_prior <- bernoulli(prob = rep(model_prior$prob, d))
      } else if (length(model_prior$prob) != d) {
        cli::cli_abort("For {.arg model_prior} of class {.cls bernoulli}, {.arg prob} must be a single value or a vector of length {.arg d}.")
      }
    }

    if (is.null(parameter_prior))
      cli::cli_abort("Argument {.arg parameter_prior} must be provided when {.arg sticky} is {.code TRUE}. It must be a single value or a vector of length {.arg d}.")

    validate_type(parameter_prior, type = "double", n = d, positive = TRUE)

  }

  grid_n <- cast_integer(grid_n, n = 1)
  validate_type(grid_n, type = "integer", n = 1, positive = TRUE)
  validate_type(grid_t_max, type = "double", n = 1, positive = TRUE)
  validate_type(post_warmup_simplify, type = "logical", n = 1)

  return(list(
    d = d, flow = flow, algorithm = algorithm, T = T, t0 = t0, t_warmup = t_warmup,
    flow_mean = flow_mean, flow_cov = flow_cov, c0 = c0,
    x0 = x0, theta0 = theta0, show_progress = show_progress,
    sticky = sticky, can_stick = can_stick, model_prior = model_prior, parameter_prior = parameter_prior,
    grid_n = grid_n, grid_t_max = grid_t_max,
    post_warmup_simplify = post_warmup_simplify,
    n_chains = n_chains, threaded = threaded, seed = seed,
    adaptive_scheme = adaptive_scheme
  ))
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
#'   during warmup and requires `"GridThinningStrategy"` as the algorithm.
#'   The `"PreconditionedZigZag"` and `"PreconditionedBPS"` flows learn a
#'   diagonal preconditioner during warmup and also require
#'   `"GridThinningStrategy"` as the algorithm.
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
#' @param sticky Logical, whether to use sticky sampling (default: FALSE).
#' @param can_stick Logical vector of length d, which coordinates can stick (default: all FALSE).
#' @param model_prior Prior distribution object for model selection. Should be of class
#'   'bernoulli' or 'beta-bernoulli' (default: NULL).
#' @param parameter_prior Numeric vector of length d, prior parameters for sticky sampling (default: NULL).
#' @param grid_n Integer, number of grid points for GridThinningStrategy (default: 30).
#' @param grid_t_max Numeric, maximum time for grid in GridThinningStrategy (default: 2.0).
#' @param post_warmup_simplify Logical. If \code{TRUE}, the grid-thinning
#'   strategy may switch to a constant bound after warmup, reducing
#'   gradient calls during the main sampling phase.
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
                        algorithm = c("ThinningStrategy", "GridThinningStrategy", "RootsPoissonStrategy"),
                        T = 50000, t0 = 0.0, t_warmup = 0.0,
                        flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                        x0 = NULL, theta0 = NULL,
                        hessian = NULL,
                        sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                        grid_n = 30, grid_t_max = 2.0,
                        post_warmup_simplify = FALSE,
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
                                sticky, can_stick, model_prior, parameter_prior,
                                grid_n, grid_t_max, post_warmup_simplify,
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

  result <- .pdmpsamplers_julia_eval("r_pdmp_custom(
    grad!, d, x0, flow, algorithm, flow_mean, flow_cov;
    c0 = c0, grid_n = grid_n, grid_t_max = grid_t_max,
    post_warmup_simplify = post_warmup_simplify,
    t0 = t0, T = T, t_warmup = t_warmup,
    hessian = hessian_f,
    sticky = sticky, can_stick = can_stick,
    model_prior = model_prior, parameter_prior = parameter_prior,
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
#' @inheritParams pdmp_sample
#'
#' @return A \code{pdmp_result} object. Use \code{mean}, \code{var},
#'   \code{quantile}, etc. for continuous-time estimators, or
#'   \code{discretize} to obtain a sample matrix.
#'
#' @export
pdmp_sample_from_stanmodel <- function(path_to_stanmodel, standata,
                        flow = c("ZigZag", "BouncyParticle", "Boomerang", "AdaptiveBoomerang", "PreconditionedZigZag", "PreconditionedBPS"),
                        algorithm = c("ThinningStrategy", "GridThinningStrategy", "RootsPoissonStrategy"),
                        T = 50000, t0 = 0.0, t_warmup = 0.0,
                        flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                        x0 = NULL, theta0 = NULL,
                        sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                        grid_n = 30, grid_t_max = 2.0,
                        post_warmup_simplify = FALSE,
                        show_progress = TRUE,
                        n_chains = 1L, threaded = FALSE, seed = NULL,
                        adaptive_scheme = c("diagonal", "fullrank"),
                        materialize = TRUE,
                        support_boundary = support_boundary_control()) {

  # Validate file paths on R side before setting up Julia
  validate_type(path_to_stanmodel, type = "character", n = 1)

  if (is.list(standata)) {
    standata_path <- tempfile(fileext = ".json")
    write_stan_json(standata, standata_path)
    on.exit(unlink(standata_path), add = TRUE)
  } else {
    validate_type(standata, type = "character", n = 1)
    standata_path <- standata
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

  check_for_julia_setup()

  # Normalize paths to absolute
  path_to_stanmodel <- normalizePath(path_to_stanmodel, mustWork = TRUE)
  standata_path     <- normalizePath(standata_path,     mustWork = TRUE)

  # Create PDMPModel in Julia and get dimension
  JuliaCall::julia_assign("_path_to_stan_model", path_to_stanmodel)
  JuliaCall::julia_assign("_path_to_stan_data",  standata_path)
  JuliaCall::julia_command("_pdmp_model = PDMPModel(_path_to_stan_model, _path_to_stan_data);")
  d <- JuliaCall::julia_eval("_pdmp_model.d")

  # Use common validation function
  params <- validate_pdmp_params(d, flow, algorithm, T, t0, t_warmup, flow_mean, flow_cov,
                                 c0, x0, theta0, show_progress,
                                 sticky, can_stick, model_prior, parameter_prior,
                                 grid_n, grid_t_max, post_warmup_simplify,
                                 n_chains, threaded, seed,
                                 adaptive_scheme = adaptive_scheme)

  support_boundary <- validate_support_boundary_control(support_boundary)

  # Pass arguments to Julia
  for (nm in names(params))
    JuliaCall::julia_assign(nm, params[[nm]])
  JuliaCall::julia_assign("support_boundary_mode", support_boundary$mode)
  JuliaCall::julia_assign("support_boundary_max_bisection_steps", support_boundary$max_bisection_steps)
  JuliaCall::julia_assign("support_boundary_time_rtol", support_boundary$time_rtol)
  JuliaCall::julia_assign("support_boundary_time_atol", support_boundary$time_atol)
  JuliaCall::julia_assign("support_boundary_clip_fraction", support_boundary$clip_fraction)
  JuliaCall::julia_assign("support_boundary_max_refresh_attempts", support_boundary$max_refresh_attempts)
  JuliaCall::julia_assign("support_boundary_refresh_probe_time", support_boundary$refresh_probe_time)
  JuliaCall::julia_assign("support_boundary_min_safe_time", support_boundary$min_safe_time)

  result <- .pdmpsamplers_julia_eval("r_pdmp_stan(
    _pdmp_model, x0, flow, algorithm, flow_mean, flow_cov;
    c0 = c0, grid_n = grid_n, grid_t_max = grid_t_max,
    post_warmup_simplify = post_warmup_simplify,
    t0 = t0, T = T, t_warmup = t_warmup,
    sticky = sticky, can_stick = can_stick,
    model_prior = model_prior, parameter_prior = parameter_prior,
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
