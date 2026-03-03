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
#' @param flow_mean Numeric vector for the flow reference mean, or NULL.
#' @param flow_cov Numeric matrix for the flow covariance, or NULL.
#' @param c0 Numeric thinning bound constant.
#' @param grid_n Integer number of grid points for GridThinningStrategy.
#' @param grid_t_max Numeric max grid interval for GridThinningStrategy.
#' @param show_progress Logical; show sampling progress bar.
#' @param discretize_dt Numeric time step for discretization, or NULL
#'   for automatic (yields ~1000 samples).
#' @param stanvars Optional `stanvar` object for custom Stan code.
#' @param sample_prior Currently only `"no"` is supported.
#' @param save_model Optional file path to save the generated Stan code.
#' @param ... Additional arguments passed to [brms::brm()] for model setup.
#'
#' @return A `brmsfit` object.
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
#' With a single chain, `Rhat` will report `NA`.
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

  empty_fit <- brms::brm(formula, data = data, family = family,
                         prior = prior, stanvars = stanvars,
                         sample_prior = sample_prior,
                         empty = TRUE, ...)

  stan_file <- tempfile(fileext = ".stan")
  cat(scode, file = stan_file)
  data_file <- tempfile(fileext = ".json")
  write_stan_json(sdata, data_file)

  if (!is.null(save_model))
    cat(scode, file = save_model)

  csv_file <- tempfile(fileext = ".csv")

  jl_flow_mean <- if (is.null(flow_mean)) numeric(0) else flow_mean
  jl_flow_cov  <- if (is.null(flow_cov)) matrix(numeric(0), nrow = 0, ncol = 0) else flow_cov
  jl_discretize_dt <- if (is.null(discretize_dt)) 0.0 else discretize_dt

  csv_paths <- JuliaCall::julia_call(
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
    n_chains = 1L,
    threaded = FALSE
  )

  stanfit <- rstan::read_stan_csv(csv_paths)

  empty_fit$fit <- stanfit
  empty_fit <- brms::rename_pars(empty_fit)

  empty_fit
}
