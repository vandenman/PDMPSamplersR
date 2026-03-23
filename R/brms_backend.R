#' Fit a brms model using PDMP samplers
#'
#' Uses PDMPSamplers.jl as a sampling backend for brms models.
#' Returns a standard `brmsfit` object with all post-processing
#' (`summary`, `plot`, `conditional_effects`, `loo`, etc.) working.
#'
#' @param formula A brms model formula.
#' @param data A data frame containing the variables in the model.
#' @param family A family object (e.g., `gaussian()`, `bernoulli()`).
#' @param prior A `brmsprior` object or NULL for default priors.
#' @param flow Character string specifying the PDMP flow type.
#' @param algorithm Character string specifying the Poisson time strategy.
#' @param adaptive_scheme Character string for preconditioner adaptation.
#' @param T Numeric total simulation time.
#' @param t0 Numeric start time.
#' @param t_warmup Numeric warmup duration. Auto-set for adaptive flows.
#'    When subsampling is active and `t_warmup` is 0, it is automatically
#'   set to 20\% of the sampling time.
#' @param flow_mean Numeric vector for the flow reference mean, or NULL.
#' @param flow_cov Numeric matrix for the flow covariance, or NULL.
#' @param c0 Numeric thinning bound constant.
#' @param grid_n Integer number of grid points for GridThinningStrategy.
#' @param grid_t_max Numeric max grid interval for GridThinningStrategy.
#' @param show_progress Logical; show sampling progress bar.
#' @param discretize_dt Numeric time step for discretization, or NULL
#'   for automatic (yields ~1000 samples).
#' @param n_chains Integer number of chains (default: 1). With multiple
#'   chains, `Rhat` and multi-chain diagnostics become available.
#' @param threaded Logical; run chains in parallel (default: FALSE).
#' @param compute_lp Logical; compute `lp__` via `BridgeStan::log_density()`
#'   for each sample (default: FALSE). Adds overhead but enables
#'   `bridge_sampler()` and populates the `lp__` diagnostic column.
#' @param subsample_size Integer number of observations per subsample,
#'   or NULL (default) for full-data gradients. When non-NULL, a
#'   BridgeStan control-variate subsampled gradient is used. Must be
#'   less than `nrow(data)`. Currently only fixed-effects models are
#'   supported for subsampling.
#' @param n_anchor_updates Integer number of anchor updates during warmup
#'   (default: 10). Only used when `subsample_size` is non-NULL.
#' @param resample_dt Numeric time step for resampling, or NULL (default)
#'   for no resampling. Only used when `subsample_size` is non-NULL.
#' @param hvp_mode Character string controlling Hessian-vector product
#'   scaling in subsampled gradients. One of `"scaled"` (default) or
#'   `"none"`.
#' @param use_hcv Logical; enable damped Hessian control variate (HCV)
#'   correction for subsampled gradients (default: FALSE). Requires
#'   `subsample_size` to be set. Adds a second-order correction that
#'   reduces gradient variance near the anchor.
#' @param use_anchor_bank Logical; enable anchor bank with LRU eviction
#'   for subsampled gradients (default: FALSE). Requires `subsample_size`
#'   to be set. Maintains multiple cached anchor points and selects the
#'   nearest one at each trajectory step.
#' @param bank_capacity Integer capacity of the anchor bank (default: 20).
#'   Only used when `use_anchor_bank` is TRUE.
#' @param use_fd_hvp Logical; use finite-difference directional curvature
#'   instead of BridgeStan's Hessian-vector product (default: FALSE).
#'   Replaces each HVP call (2-3x gradient cost) with one extra gradient
#'   call, giving ~30-50\% per-grid-point savings.
#' @param post_warmup_simplify Logical; switch to a fast constant-bound
#'   thinning mode after warmup when the sampler is well-adapted
#'   (default: FALSE). Activates when the reflection ratio is below 30\%
#'   and at least 10 events have been observed.
#' @param use_fd_hcv Logical; use finite-difference approximation for the
#'   HVP inside the HCV correction (default: FALSE). Replaces the exact
#'   BridgeStan HVP at the anchor with a cheaper gradient-based FD
#'   approximation. Also disables HVP for grid thinning. Only used when
#'   `use_hcv` is TRUE.
#' @param stanvars Optional `stanvar` object for custom Stan code.
#' @param sample_prior Currently only `"no"` is supported.
#' @param save_model Optional file path to save the generated Stan code.
#' @param ... Additional arguments passed to [brms::brm()] for model setup.
#'
#' @return A `brmsfit` object.
#'
#' @importFrom stats gaussian
#'
#' @details
#' PDMP samplers do not produce NUTS-style diagnostics. The `lp__`,
#' `accept_stat__`, and related diagnostic columns are set to zero.
#' Functions relying on NUTS diagnostics (e.g., `pairs()` divergence
#' plots, `nuts_params()`) will not produce meaningful output.
#'
#' `loo()` works because brms computes `log_lik` directly from model
#' parameters, not from `lp__`.
#'
#' With a single chain, `Rhat` will report `NA`. Use `n_chains >= 2`
#' for convergence diagnostics.
#'
#' When `subsample_size` is specified, the function uses BridgeStan
#' data-swapping to compute control-variate subsampled gradients. A
#' centering fix is applied to the brms-generated Stan code so that
#' predictor centering remains consistent across subsamples.
#'
#' @export
brm_pdmp <- function(
    formula, data, family = gaussian(),
    prior = NULL,
    flow = c("ZigZag", "BouncyParticle", "Boomerang",
             "AdaptiveBoomerang", "PreconditionedZigZag", "PreconditionedBPS"),
    algorithm = c("GridThinningStrategy", "ThinningStrategy",
                  "RootsPoissonStrategy"),
    adaptive_scheme = c("diagonal", "fullrank"),
    T = 50000, t0 = 0.0, t_warmup = 0.0,
    flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
    grid_n = 30, grid_t_max = 2.0,
    show_progress = TRUE,
    discretize_dt = NULL,
    n_chains = 1L, threaded = FALSE,
    compute_lp = FALSE,
    subsample_size = NULL,
    n_anchor_updates = 10L,
    resample_dt = NULL,
    hvp_mode = c("scaled", "none"),
    use_hcv = FALSE,
    use_anchor_bank = FALSE,
    bank_capacity = 20L,
    use_fd_hvp = FALSE,
    post_warmup_simplify = FALSE,
    use_fd_hcv = FALSE,
    stanvars = NULL, sample_prior = "no",
    save_model = NULL,
    ...
) {
  if (!requireNamespace("brms", quietly = TRUE))
    cli::cli_abort("Package {.pkg brms} is required for {.fn brm_pdmp}.")
  if (!requireNamespace("rstan", quietly = TRUE))
    cli::cli_abort("Package {.pkg rstan} is required for {.fn brm_pdmp}.")

  if (sample_prior != "no")
    cli::cli_abort("{.fn brm_pdmp} only supports {.code sample_prior = \"no\"}.")

  flow <- match.arg(flow)
  algorithm <- match.arg(algorithm)
  adaptive_scheme <- match.arg(adaptive_scheme)
  hvp_mode <- match.arg(hvp_mode)
  subsampled <- !is.null(subsample_size)
  N <- nrow(data)

  if (!subsampled && (use_hcv || use_anchor_bank))
    cli::cli_abort("{.arg use_hcv} and {.arg use_anchor_bank} require {.arg subsample_size} to be set.")

  if (subsampled) {
    bank_capacity <- as.integer(bank_capacity)
    subsample_size <- as.integer(subsample_size)
    if (subsample_size >= N)
      cli::cli_abort("{.arg subsample_size} ({subsample_size}) must be less than {.code nrow(data)} ({N}).")
    n_anchor_updates <- as.integer(n_anchor_updates)
    if (t_warmup == 0) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for subsampled gradients.")
    }
  }

  if (flow == "AdaptiveBoomerang") {
    if (algorithm != "GridThinningStrategy")
      cli::cli_abort("{.val AdaptiveBoomerang} requires {.val GridThinningStrategy} as the algorithm.")
    if (t_warmup == 0) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for {.val AdaptiveBoomerang}.")
    }
  }
  if (flow %in% c("PreconditionedZigZag", "PreconditionedBPS")) {
    if (algorithm != "GridThinningStrategy")
      cli::cli_abort("{.val {flow}} requires {.val GridThinningStrategy} as the algorithm.")
    if (t_warmup == 0) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for {.val {flow}}.")
    }
  }

  check_for_julia_setup()

  scode <- brms::stancode(formula, data = data, family = family,
                          prior = prior, stanvars = stanvars,
                          sample_prior = sample_prior, ...)
  sdata <- brms::standata(formula, data = data, family = family,
                          prior = prior, stanvars = stanvars,
                          sample_prior = sample_prior, ...)

  if (subsampled) {
    if (any(grepl("^(J|Z)_", names(sdata))))
      cli::cli_abort(c(
        "Subsampled gradients currently support fixed-effects models only.",
        "i" = "Remove random effects from the formula, or omit {.arg subsample_size} to use full-data gradients."
      ))
    scode <- fix_brms_stancode(scode)
    scode_ext <- inject_ext_cpp_stancode(scode)
    if (!is.null(sdata$means_X)) {
      means_X <- sdata$means_X
    } else {
      means_X <- array(colMeans(sdata$X[, -1, drop = FALSE]))
    }
    Y_full <- as.numeric(sdata$Y)
    X_full <- sdata$X
    sdata_full <- sdata
    sdata_full$means_X <- means_X
    sdata_prior <- make_prior_standata(sdata, means_X)
  }

  empty_fit <- brms::brm(formula, data = data, family = family,
                         prior = prior, stanvars = stanvars,
                         sample_prior = sample_prior,
                         empty = TRUE, ...)

  stan_file <- cached_stan_model(scode)

  if (subsampled) {
    stan_file_ext <- cached_stan_model(scode_ext)
    data_full_file <- tempfile(fileext = ".json")
    write_stan_json(sdata_full, data_full_file)
    data_prior_file <- tempfile(fileext = ".json")
    write_stan_json(sdata_prior, data_prior_file)
  } else {
    data_file <- tempfile(fileext = ".json")
    write_stan_json(sdata, data_file)
  }

  if (!is.null(save_model))
    cat(scode, file = save_model)

  csv_file <- tempfile(fileext = ".csv")
  jl_flow_mean <- if (is.null(flow_mean)) numeric(0) else flow_mean
  jl_flow_cov  <- if (is.null(flow_cov)) matrix(numeric(0), nrow = 0, ncol = 0) else flow_cov
  jl_discretize_dt <- if (is.null(discretize_dt)) 0.0 else discretize_dt
  jl_resample_dt <- if (is.null(resample_dt)) 0.0 else resample_dt

  if (subsampled) {
    jl_result <- JuliaCall::julia_call(
      "r_pdmp_brms_subsampled",
      normalizePath(stan_file, mustWork = TRUE),
      normalizePath(stan_file_ext, mustWork = TRUE),
      normalizePath(hpp_path(), mustWork = TRUE),
      normalizePath(data_full_file, mustWork = TRUE),
      normalizePath(data_prior_file, mustWork = TRUE),
      as.integer(N), subsample_size,
      flow, algorithm,
      jl_flow_mean, jl_flow_cov,
      csv_file,
      c0 = c0,
      grid_n = as.integer(grid_n),
      grid_t_max = grid_t_max,
      t0 = t0, T = T, t_warmup = t_warmup,
      n_anchor_updates = n_anchor_updates,
      adaptive_scheme = adaptive_scheme,
      discretize_dt = jl_discretize_dt,
      show_progress = show_progress,
      n_chains = as.integer(n_chains),
      threaded = threaded,
      compute_lp = compute_lp,
      resample_dt = jl_resample_dt,
      hvp_mode = hvp_mode,
      use_hcv = use_hcv,
      use_anchor_bank = use_anchor_bank,
      bank_capacity = as.integer(bank_capacity),
      use_fd_hvp = use_fd_hvp,
      post_warmup_simplify = post_warmup_simplify,
      use_fd_hcv = use_fd_hcv
    )
  } else {
    jl_result <- JuliaCall::julia_call(
      "r_pdmp_stan_for_brms",
      normalizePath(stan_file, mustWork = TRUE),
      normalizePath(data_file, mustWork = TRUE),
      flow, algorithm,
      jl_flow_mean, jl_flow_cov,
      csv_file,
      c0 = c0,
      grid_n = as.integer(grid_n),
      grid_t_max = grid_t_max,
      t0 = t0, T = T, t_warmup = t_warmup,
      adaptive_scheme = adaptive_scheme,
      discretize_dt = jl_discretize_dt,
      show_progress = show_progress,
      n_chains = as.integer(n_chains),
      threaded = threaded,
      compute_lp = compute_lp,
      use_fd_hvp = use_fd_hvp,
      post_warmup_simplify = post_warmup_simplify
    )
  }

  if (is.environment(jl_result)) jl_result <- as.list(jl_result)
  csv_paths <- jl_result$csv_paths
  pdmp_stats <- jl_result$stats
  sampling_time <- pdmp_stats$elapsed_time

  stanfit <- rstan::read_stan_csv(csv_paths)
  empty_fit$fit <- stanfit
  empty_fit <- brms::rename_pars(empty_fit)
  attr(empty_fit, "sampling_time") <- sampling_time
  attr(empty_fit, "pdmp_stats") <- pdmp_stats
  empty_fit
}

cached_stan_model <- function(scode) {
  cache_dir <- file.path(rappdirs::user_cache_dir("PDMPSamplersR"), "stan_cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  hash <- rlang::hash(scode)
  model_base <- file.path(cache_dir, paste0(hash, "_model"))
  lib_candidates <- c(
    paste0(model_base, ".so"),
    paste0(model_base, ".dylib"),
    paste0(model_base, ".dll")
  )
  existing_lib <- lib_candidates[file.exists(lib_candidates)]
  if (length(existing_lib) > 0) return(existing_lib[[1]])
  stan_path <- file.path(cache_dir, paste0(hash, ".stan"))
  writeLines(scode, stan_path)
  stan_path
}
