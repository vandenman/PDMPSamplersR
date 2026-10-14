#' Additional sampler controls
#'
#' Pass this object as `control` to [pdmp_sample()],
#' [pdmp_sample_from_stanmodel()] or [brm_pdmp()]. Settings apply to that call
#' only, including when a worker is reused for another model.
#'
#' @param grid_n_min Minimum grid size, or `NULL` for `min(grid_n, 5)`.
#' @param grid_bound_violation Grid bound violation policy: `NULL` uses the
#'   sampler default; otherwise `"throw"`, `"shrink"`, or `"count"`.
#' @param adaptive_boomerang_refresh_rate Initial refresh rate for adaptive
#'   Boomerang dynamics. The default is 0.1.
#' @param warmup_stop Optional named list configuring adaptive warmup stopping:
#'   `min_time`, `max_time`, `stable_time`, `min_events`, and `check_every`.
#'   An empty list enables defaults (one quarter, all, and one quarter of
#'   `t_warmup`, 100 events, and a check every 25 events); `NULL` disables it.
#'   Supported by custom-Stan sampling.
#' @param profile Collect bridge-phase timings for custom-Stan sampling.
#' @param research_diagnostics Collect additional OMRF allocation and HVP timing
#'   diagnostics. These measurements add overhead.
#' @param hcv_cache Reuse OMRF Hessian products along a flight.
#' @param check_fast_deterministic Compare the fast OMRF deterministic gradient
#'   with its generic evaluation. Intended for diagnostics; adds overhead.
#' @return A sampler control list.
#' @md
#' @export
sampler_control <- function(grid_n_min = NULL, grid_bound_violation = NULL,
    adaptive_boomerang_refresh_rate = 0.1, warmup_stop = NULL,
    profile = FALSE, research_diagnostics = FALSE, hcv_cache = TRUE,
    check_fast_deterministic = FALSE) {
  if (!is.null(grid_n_min) &&
      (!rlang::is_integerish(grid_n_min, n = 1L, finite = TRUE) || grid_n_min < 1))
    cli::cli_abort("{.arg grid_n_min} must be a positive integer or NULL.")
  if (!is.null(grid_bound_violation))
    grid_bound_violation <- match.arg(grid_bound_violation, c("throw", "shrink", "count"))
  if (!is.numeric(adaptive_boomerang_refresh_rate) ||
      length(adaptive_boomerang_refresh_rate) != 1L ||
      !is.finite(adaptive_boomerang_refresh_rate) || adaptive_boomerang_refresh_rate < 0)
    cli::cli_abort("{.arg adaptive_boomerang_refresh_rate} must be finite and nonnegative.")
  for (name in c("profile", "research_diagnostics", "hcv_cache", "check_fast_deterministic")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      cli::cli_abort("{.arg {name}} must be TRUE or FALSE.")
  }
  if (!is.null(warmup_stop)) {
    allowed <- c("min_time", "max_time", "stable_time", "min_events", "check_every")
    if (!is.list(warmup_stop) || (length(warmup_stop) &&
        (is.null(names(warmup_stop)) || anyDuplicated(names(warmup_stop)) ||
         any(!names(warmup_stop) %in% allowed))))
      cli::cli_abort("{.arg warmup_stop} must be NULL or a named list of warmup stopping settings.")
    for (name in names(warmup_stop)) {
      value <- warmup_stop[[name]]
      integer <- name %in% c("min_events", "check_every")
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
          value < 0 || (integer && (!rlang::is_integerish(value) || value < 1)))
        cli::cli_abort("Invalid warmup stopping setting {.arg {name}}.")
    }
  }
  structure(list(grid_n_min = if (is.null(grid_n_min)) NULL else as.integer(grid_n_min),
    grid_bound_violation = grid_bound_violation,
    adaptive_boomerang_refresh_rate = as.numeric(adaptive_boomerang_refresh_rate),
    warmup_stop = warmup_stop, profile = profile,
    research_diagnostics = research_diagnostics, hcv_cache = hcv_cache,
    check_fast_deterministic = check_fast_deterministic), class = "pdmp_sampler_control")
}

.validate_sampler_control <- function(control) {
  if (!inherits(control, "pdmp_sampler_control"))
    cli::cli_abort("{.arg control} must come from {.fn sampler_control}.")
  do.call(sampler_control, unclass(control))
}

#' Stan compilation options
#'
#' Pass as `compile_options` to [pdmp_sample_from_stanmodel()], [brm_pdmp()]
#' or [compile_pdmp_stan_model()]. Existing compiled libraries can be reused;
#' use a separate cache directory when changing toolchains.
#' @param cache Use a content-addressed compilation cache.
#' @param cache_dir Cache directory, or `NULL` for the Julia depot's
#'   `pdmpsamplersr/bridgestan` directory.
#' @param native Enable native CPU compiler instructions.
#' @param optimize Enable Stan's `--O1` optimization.
#' @return A compilation options list.
#' @md
#' @export
stan_compile_options <- function(cache = FALSE, cache_dir = NULL,
                                 native = FALSE, optimize = FALSE) {
  for (name in c("cache", "native", "optimize")) {
    value <- get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value))
      cli::cli_abort("{.arg {name}} must be TRUE or FALSE.")
  }
  if (!is.null(cache_dir)) {
    validate_type(cache_dir, type = "character", n = 1)
    if (is.na(cache_dir) || !nzchar(cache_dir)) cli::cli_abort("{.arg cache_dir} must not be empty.")
    cache_dir <- normalizePath(path.expand(cache_dir), mustWork = FALSE)
  }
  structure(list(cache = cache, cache_dir = cache_dir, native = native,
    optimize = optimize), class = "pdmp_stan_compile_options")
}

.validate_compile_options <- function(compile_options) {
  if (!inherits(compile_options, "pdmp_stan_compile_options"))
    cli::cli_abort("{.arg compile_options} must come from {.fn stan_compile_options}.")
  do.call(stan_compile_options, unclass(compile_options))
}
