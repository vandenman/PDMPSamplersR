#' PDMP sampling with subsampled gradients
#'
#' Performs PDMP sampling using a user-supplied subsampled gradient callback.
#' This is a low-level interface intended for correctness verification and
#' prototyping. For large-scale problems with cheap likelihoods, consider the
#' Julia-native built-in targets (when available) to avoid R-to-Julia call
#' overhead.
#'
#' @param grad_sub Function that computes the NEGATIVE subsampled gradient.
#'   Should take two arguments: a numeric vector `x` (position) and an integer
#'   vector `indices` (1-based observation indices), and return a numeric vector
#'   of length `d`.
#' @param n_obs Integer, total number of observations in the dataset.
#' @param d Integer, dimension of the parameter space.
#' @param subsample_size Integer, number of observations per subsample.
#' @param flow Character string specifying the flow type. One of `"ZigZag"`,
#'   `"BouncyParticle"`, `"Boomerang"`, `"AdaptiveBoomerang"`,
#'   `"PreconditionedZigZag"`, or `"PreconditionedBPS"`.
#' @param algorithm Character string specifying the algorithm. One of
#'   `"ThinningStrategy"`, `"GridThinningStrategy"`, or
#'   `"RootsPoissonStrategy"`.
#' @param T Numeric, total sampling time (default: 50000).
#' @param t0 Numeric, initial time (default: 0.0).
#' @param t_warmup Numeric, warmup time (default: 0.0).
#' @param flow_mean Numeric vector of length d, mean for the flow (default: zero vector).
#' @param flow_cov Numeric matrix of size d x d, covariance for the flow (default: identity).
#' @param c0 Numeric, bound parameter (default: 1e-2).
#' @param x0 Numeric vector of length d, initial position (default: random normal).
#' @param hvp_sub Function for the subsampled Hessian-vector product (optional).
#'   Takes `(x, v, indices)` and returns a numeric vector of length `d`.
#'   If `NULL`, a finite-difference approximation is used.
#' @param grad_full Function for the full (non-subsampled) NEGATIVE gradient (optional).
#'   Takes a numeric vector `x` and returns a numeric vector of length `d`.
#'   Required when `use_full_gradient_for_reflections = TRUE`.
#' @param use_full_gradient_for_reflections Logical, whether to use the full gradient
#'   for reflection events (default: `FALSE`).
#' @param grid_n Integer, number of grid points for GridThinningStrategy (default: 30).
#' @param grid_t_max Numeric, maximum time for grid (default: 2.0).
#' @param show_progress Logical, whether to show progress bar (default: `TRUE`).
#' @param n_chains Integer, number of chains to run (default: 1).
#' @param threaded Logical, whether to run chains in parallel (default: `FALSE`).
#' @param adaptive_scheme Character string, adaptation scheme for AdaptiveBoomerang
#'   (default: `"diagonal"`).
#'
#' @return A \code{pdmp_result} object.
#'
#' @export
pdmp_sample_subsampled <- function(
    grad_sub, n_obs, d, subsample_size,
    flow = c("ZigZag", "BouncyParticle", "Boomerang",
             "AdaptiveBoomerang", "PreconditionedZigZag", "PreconditionedBPS"),
    algorithm = c("GridThinningStrategy", "ThinningStrategy",
                  "RootsPoissonStrategy"),
    T = 50000, t0 = 0.0, t_warmup = 0.0,
    flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
    x0 = NULL,
    hvp_sub = NULL,
    grad_full = NULL,
    use_full_gradient_for_reflections = FALSE,
    grid_n = 30, grid_t_max = 2.0,
    show_progress = TRUE,
    n_chains = 1L, threaded = FALSE,
    adaptive_scheme = c("diagonal", "fullrank")) {

  if (!rlang::is_function(grad_sub))
    cli::cli_abort("{.arg grad_sub} must be a function, not {.cls {class(grad_sub)}}.")

  n_obs <- cast_integer(n_obs, n = 1)
  validate_type(n_obs, type = "integer", n = 1, positive = TRUE)

  subsample_size <- cast_integer(subsample_size, n = 1)
  validate_type(subsample_size, type = "integer", n = 1, positive = TRUE)

  if (subsample_size >= n_obs)
    cli::cli_abort("{.arg subsample_size} ({subsample_size}) must be less than {.arg n_obs} ({n_obs}).")

  params <- validate_pdmp_params(
    d, flow, algorithm, T, t0, t_warmup, flow_mean, flow_cov,
    c0, x0, theta0 = NULL, show_progress,
    sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
    grid_n, grid_t_max, n_chains, threaded,
    adaptive_scheme = adaptive_scheme
  )

  test_indices <- sample.int(n_obs, subsample_size)
  tryCatch({
    test_out <- grad_sub(params$x0, test_indices)
    if (!is.numeric(test_out) || length(test_out) != params$d)
      cli::cli_abort("{.arg grad_sub} must return a numeric vector of length {params$d}, got length {length(test_out)}.")
  }, error = function(e) {
    cli::cli_abort(c(
      "{.arg grad_sub} failed on test input.",
      "x" = conditionMessage(e)
    ))
  })

  if (!is.null(hvp_sub)) {
    if (!rlang::is_function(hvp_sub))
      cli::cli_abort("{.arg hvp_sub} must be a function or NULL, not {.cls {class(hvp_sub)}}.")
    tryCatch({
      test_v <- stats::rnorm(params$d)
      test_out <- hvp_sub(params$x0, test_v, test_indices)
      if (!is.numeric(test_out) || length(test_out) != params$d)
        cli::cli_abort("{.arg hvp_sub} must return a numeric vector of length {params$d}, got length {length(test_out)}.")
    }, error = function(e) {
      cli::cli_abort(c(
        "{.arg hvp_sub} failed on test input.",
        "x" = conditionMessage(e)
      ))
    })
  }

  if (!is.null(grad_full)) {
    if (!rlang::is_function(grad_full))
      cli::cli_abort("{.arg grad_full} must be a function or NULL, not {.cls {class(grad_full)}}.")
    tryCatch({
      test_out <- grad_full(params$x0)
      if (!is.numeric(test_out) || length(test_out) != params$d)
        cli::cli_abort("{.arg grad_full} must return a numeric vector of length {params$d}, got length {length(test_out)}.")
    }, error = function(e) {
      cli::cli_abort(c(
        "{.arg grad_full} failed on test input.",
        "x" = conditionMessage(e)
      ))
    })
  }

  if (use_full_gradient_for_reflections && is.null(grad_full))
    cli::cli_abort("{.arg grad_full} must be provided when {.arg use_full_gradient_for_reflections} is {.code TRUE}.")

  check_for_julia_setup()

  for (nm in names(params))
    JuliaCall::julia_assign(nm, params[[nm]])

  JuliaCall::julia_assign("_n_obs", n_obs)
  JuliaCall::julia_assign("_subsample_size", subsample_size)
  JuliaCall::julia_assign("_grad_sub_r", grad_sub)
  JuliaCall::julia_assign("_hvp_sub_r", hvp_sub)
  JuliaCall::julia_assign("_grad_full_r", grad_full)
  JuliaCall::julia_assign("_use_full_for_reflections", use_full_gradient_for_reflections)

  result <- JuliaCall::julia_eval("r_pdmp_custom_subsampled(
    _grad_sub_r, d, _n_obs, _subsample_size, x0, flow, algorithm,
    flow_mean, flow_cov;
    hvp_sub_r = _hvp_sub_r,
    grad_full_r = _grad_full_r,
    use_full_gradient_for_reflections = _use_full_for_reflections,
    c0 = c0, grid_n = grid_n, grid_t_max = grid_t_max,
    t0 = t0, T = T, t_warmup = t_warmup,
    show_progress = show_progress, n_chains = n_chains, threaded = threaded,
    adaptive_scheme = adaptive_scheme
  );")
  if (is.environment(result)) result <- as.list(result)
  new_pdmp_result(
    chains   = result$chains,
    stats    = result$stats,
    d        = result$d,
    n_chains = result$n_chains
  )
}
