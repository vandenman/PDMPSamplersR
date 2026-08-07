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
#' @param adaptive_scheme Character string for AdaptiveBoomerang covariance
#'   adaptation: `"diagonal"`, `"fullrank"`, or `"lowrank"`.
#' @param T Numeric total simulation time.
#' @param t0 Numeric start time.
#' @param t_warmup Numeric warmup duration. Auto-set for adaptive flows.
#'    When subsampling is active and `t_warmup` is 0, it is automatically
#'   set to 20\% of the sampling time.
#' @param flow_mean Numeric vector for the flow reference mean, or NULL. For
#'   marked subsampling this is also the fixed control-variate anchor and
#'   initial position; supplying a posterior mode or other representative point
#'   can materially tighten the residual envelope.
#' @param flow_cov Numeric matrix for the flow covariance, or NULL.
#' @param c0 Numeric thinning bound constant. With `ThinningStrategy` this
#'   user-configured bound must dominate the deterministic marked component;
#'   violations are detected and stop the sampler.
#' @param grid_n Integer number of grid points for GridThinningStrategy.
#' @param grid_t_max Numeric max grid interval for GridThinningStrategy.
#' @param show_progress Logical; show sampling progress bar.
#' @param discretize_dt Numeric time step for discretization, or NULL
#'   for automatic (yields ~1000 samples).
#' @param n_chains Integer number of chains (default: 1). With multiple
#'   chains, `Rhat` and multi-chain diagnostics become available.
#' @param threaded Logical; run chains in parallel (default: FALSE).
#' @param seed NULL (default) or a non-negative integer seed passed through to
#'   Julia's sampler RNG.
#' @param compute_lp Logical; compute `lp__` via `BridgeStan::log_density()`
#'   for each sample (default: FALSE). Adds overhead but enables
#'   `bridge_sampler()` and populates the `lp__` diagnostic column.
#' @param subsample_size Integer number of observations per subsample,
#'   or NULL (default) for full-data gradients. With `GridThinningStrategy`
#'   or a representable `ThinningStrategy` envelope,
#'   marked acceleration is selected automatically when affine predictor,
#'   first-batch likelihood, and trajectory providers are all available for
#'   the requested dynamics. Ineligible models use exact full gradients.
#' @param n_anchor_updates Integer number of anchor updates during warmup
#'   (default: 0). For analytic marked sampling, each update prepares a
#'   coherent anchor state containing the full and deterministic anchor
#'   gradients, cached likelihood predictors, and residual envelope.
#' @param use_hcv Logical; enable the certified analytic damped Hessian
#'   control variate for marked affine Bernoulli-logit or binomial-logit
#'   sampling (default: FALSE). Its likelihood Hessian and Taylor-remainder
#'   envelope are prepared with each anchor, with Cauchy damping
#'   `d / (d + ||x - anchor||^2)` in `d` unconstrained dimensions;
#'   finite-difference HCV is not used.
#' @param use_anchor_bank Logical; retain prepared warmup anchors and select
#'   the nearest one at event boundaries (default: FALSE). Selection is
#'   chain-local and performs no full-data gradient or envelope rebuild.
#'   Requires `n_anchor_updates > 0`; use only `n_anchor_updates` for the
#'   lower-memory single-anchor update mode.
#' @param bank_capacity Integer capacity of the anchor bank (default: 20).
#'   Least-recently-used entries are replaced during warmup. Each populated
#'   entry stores predictor and kernel caches of order `N * K` plus a complete
#'   residual envelope of order `N * R`; memory therefore grows linearly in
#'   the populated capacity and can be substantial for large `N` or many
#'   growth-envelope components. Benchmark representative models before using
#'   a large capacity. Analytic HCV additionally stores one dense `d * d`
#'   likelihood Hessian per populated anchor.
#' @param use_fd_hvp Logical; use finite-difference directional curvature
#'   instead of BridgeStan's Hessian-vector product (default: FALSE).
#'   Replaces each HVP call (2-3x gradient cost) with one extra gradient
#'   call, giving ~30-50\% per-grid-point savings.
#' @param post_warmup_simplify Logical; switch to a fast constant-bound
#'   thinning mode after warmup when the sampler is well-adapted
#'   (default: FALSE). Activates when the reflection ratio is below 30\%
#'   and at least 10 events have been observed.
#' @param sticky Logical; enable spike-and-slab variable selection for
#'   population-level coefficients (default: FALSE). Requires
#'   `model_prior` to be set. Sticky transitions redraw the complete velocity
#'   from the exact law on each active coordinate stratum; no pre-freeze
#'   velocity is restored. This applies to every supported dynamics, including
#'   dense-preconditioned Zig-Zag.
#' @param can_stick Optional logical vector indicating which non-intercept
#'   population-level coefficients are candidates for selection. Length
#'   must match the number of supported coefficients. If omitted, all
#'   non-intercept population-level coefficients are candidates.
#' @param model_prior A [bernoulli()] or [betabernoulli()] object specifying
#'   the prior on model space for the currently supported legacy sticky path.
#'   [exchangeable_model_size_prior()] is reserved for the pending dependent
#'   slab path.
#' @param kappa Optional numeric vector of slab densities at zero for each
#'   stickable coordinate (κ in the sticky PDMP literature). If omitted,
#'   derived automatically from the brms prior specification (only
#'   `normal(0, s)` and `student_t(df, 0, s)` are supported for automatic
#'   derivation).
#' @param slab_prior Optional dependent slab prior created by
#'   [dense_gaussian_slab()], [exchangeable_gaussian_slab()],
#'   [independent_slab_density()], [gaussian_scale_mixture_slab()], or
#'   [arbitrary_slab_boundary()]. Mutually exclusive with `kappa`. Full-data
#'   dependent slabs use the brms prior-only data as the base prior target;
#'   subsampled brms dependent slabs are not yet supported.
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
#' Marked acceleration is enabled when the model has certified affine predictor,
#' likelihood, and dynamics providers. The first provider batch covers weighted
#' and subsetted Bernoulli/binomial logit models, categorical/multinomial logit,
#' Poisson log models (including `rate()`), and Gaussian location-scale models
#' (including known `se()` values and affine sigma predictors). Independent
#' affine multivariate Bernoulli responses are also supported. Unknown geometry
#' falls back to the exact full-gradient sampler. Priors, Jacobians, and
#' unconditional custom target additions remain in an opaque deterministic
#' BridgeStan provider.
#'
#' PDMPSamplers.jl applies `N / m`, draws a fresh marked subset per proposal,
#' uses dynamics-specific linear or harmonic trajectory bounds, and reuses the
#' accepted stochastic gradient for the event. Marked subsampling supports
#' `GridThinningStrategy` generally and `ThinningStrategy` when the residual
#' envelope has an affine global roof (or the trajectory is periodic).
#'
#' @export
brm_pdmp <- function(
    formula, data, family = gaussian(),
    prior = NULL,
    flow = c("ZigZag", "BouncyParticle", "Boomerang",
             "AdaptiveBoomerang", "PreconditionedZigZag", "PreconditionedBPS",
             "DensePreconditionedZigZag", "DensePreconditionedBPS"),
    algorithm = c("GridThinningStrategy", "ThinningStrategy",
                  "RootsPoissonStrategy"),
    adaptive_scheme = c("diagonal", "fullrank", "lowrank"),
    T = 50000, t0 = 0.0, t_warmup = 0.0,
    flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
    grid_n = 30, grid_t_max = 2.0,
    show_progress = TRUE,
    discretize_dt = NULL,
    n_chains = 1L, threaded = FALSE, seed = NULL,
    compute_lp = FALSE,
    subsample_size = NULL,
    n_anchor_updates = 0L,
    use_hcv = FALSE,
    use_anchor_bank = FALSE,
    bank_capacity = 20L,
    use_fd_hvp = FALSE,
    post_warmup_simplify = FALSE,
    sticky = FALSE, can_stick = NULL, model_prior = NULL,
    kappa = NULL, slab_prior = NULL,
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
  if (!is.null(slab_prior)) {
    if (!isTRUE(sticky)) {
      cli::cli_abort("Argument {.arg slab_prior} requires {.arg sticky} to be {.code TRUE}.")
    }
    if (!flow %in% c("ZigZag", "BouncyParticle") || algorithm != "GridThinningStrategy") {
      cli::cli_abort("Dependent {.arg slab_prior} sticky sampling currently requires ZigZag or BouncyParticle with {.val GridThinningStrategy}.")
    }
    cli::cli_abort("Dependent {.arg slab_prior} is temporarily gated for {.fn brm_pdmp} until target composition subtracts only the slab component or adds back nuisance priors.")
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
  subsampled <- !is.null(subsample_size)
  requested_t_warmup <- t_warmup
  if (subsampled && !is.null(slab_prior)) {
    cli::cli_abort("Dependent {.arg slab_prior} is not yet supported together with {.arg subsample_size}.")
  }
  N <- nrow(data)

  if (!subsampled && (use_hcv || use_anchor_bank))
    cli::cli_abort("{.arg use_hcv} and {.arg use_anchor_bank} require {.arg subsample_size} to be set.")

  if (subsampled) {
    if (!rlang::is_integerish(subsample_size, n = 1L, finite = TRUE) ||
        subsample_size <= 0)
      cli::cli_abort("{.arg subsample_size} must be a positive integerish scalar.")
    if (!rlang::is_integerish(bank_capacity, n = 1L, finite = TRUE) ||
        bank_capacity <= 0)
      cli::cli_abort("{.arg bank_capacity} must be a positive integer.")
    if (!rlang::is_integerish(n_anchor_updates, n = 1L, finite = TRUE) ||
        n_anchor_updates < 0)
      cli::cli_abort("{.arg n_anchor_updates} must be a non-negative integer.")
    subsample_size <- as.integer(subsample_size)
    bank_capacity <- as.integer(bank_capacity)
    n_anchor_updates <- as.integer(n_anchor_updates)
    if (subsample_size >= N)
      cli::cli_abort("{.arg subsample_size} ({subsample_size}) must be less than {.code nrow(data)} ({N}).")
    if (isTRUE(use_anchor_bank) && n_anchor_updates == 0L)
      cli::cli_abort(paste0(
        "{.arg use_anchor_bank} requires a positive {.arg n_anchor_updates}; ",
        "otherwise no additional anchors can be prepared."
      ))
    if (t_warmup == 0) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for subsampled gradients.")
    }
  }

  if (flow == "AdaptiveBoomerang") {
    if (!algorithm %in% c("GridThinningStrategy", "ThinningStrategy"))
      cli::cli_abort("{.val AdaptiveBoomerang} requires a thinning algorithm with a certified trajectory bound.")
    if (t_warmup == 0) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for {.val AdaptiveBoomerang}.")
    }
  }
  if (flow %in% c("PreconditionedZigZag", "PreconditionedBPS",
                  "DensePreconditionedZigZag", "DensePreconditionedBPS")) {
    if (!algorithm %in% c("GridThinningStrategy", "ThinningStrategy"))
      cli::cli_abort("{.val {flow}} requires a thinning algorithm with a certified trajectory bound.")
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
    eligibility <- marked_subsampling_eligibility(formula, family, stanvars, sdata)
    runtime_supported <- flow %in% c(
      "BouncyParticle", "ZigZag", "PreconditionedBPS",
      "PreconditionedZigZag", "DensePreconditionedBPS",
      "DensePreconditionedZigZag", "Boomerang", "AdaptiveBoomerang"
    ) && algorithm %in% c("GridThinningStrategy", "ThinningStrategy")
    periodic_flow <- flow %in% c("Boomerang", "AdaptiveBoomerang")
    affine_thinning_envelope <- isTRUE(eligibility$eligible) &&
      !identical(eligibility$family, "poisson")
    if (identical(algorithm, "ThinningStrategy") && isTRUE(eligibility$eligible) &&
        !periodic_flow && !affine_thinning_envelope) {
      runtime_supported <- FALSE
      eligibility <- list(
        eligible = FALSE,
        reason = paste0(
          "ThinningStrategy has only an affine global clock, which cannot dominate ",
          "this model's exponentially growing residual envelope on a linear trajectory; ",
          "use GridThinningStrategy"
        )
      )
    }
    if (!runtime_supported) {
      if (isTRUE(eligibility$eligible)) {
        eligibility <- list(
          eligible = FALSE,
          reason = paste0(
            "the marked fast path requires a dynamics trajectory provider and ",
            "GridThinningStrategy or ThinningStrategy"
          )
        )
      }
    }
    if (!isTRUE(eligibility$eligible)) {
      cli::cli_inform(c(
        "Subsampling is not eligible for this model; using the full-gradient sampler.",
        "i" = eligibility$reason
      ))
      subsampled <- FALSE
      use_hcv <- FALSE
      use_anchor_bank <- FALSE
      if (!flow %in% c("AdaptiveBoomerang", "PreconditionedZigZag",
                       "PreconditionedBPS", "DensePreconditionedZigZag",
                       "DensePreconditionedBPS")) {
        t_warmup <- requested_t_warmup
      }
    } else {
      N <- as.integer(sdata$N)
      if (subsample_size >= N) {
        cli::cli_abort(paste0(
          "subsample_size (", subsample_size,
          ") must be smaller than the number of included observations (", N,
          ") after applying subset() modifiers."
        ))
      }
      if (isTRUE(use_hcv) && !eligibility$family %in% c("bernoulli", "binomial")) {
        cli::cli_abort(paste0(
          "The analytic marked HCV currently requires an affine Bernoulli-logit ",
          "or binomial-logit likelihood."
        ))
      }
    }
  }

  if (subsampled || !is.null(slab_prior)) {
    sdata_prior <- if (subsampled) {
      make_opaque_deterministic_standata(sdata)
    } else {
      make_prior_standata(sdata)
    }
  }

  empty_fit <- brms::brm(formula, data = data, family = family,
                         prior = prior, stanvars = stanvars,
                         sample_prior = sample_prior,
                         empty = TRUE, ...)

  stan_file <- cached_stan_model(scode)

  if (subsampled) {
    data_full_file <- tempfile(fileext = ".json")
    write_stan_json(sdata, data_full_file)
    data_prior_file <- tempfile(fileext = ".json")
    write_stan_json(sdata_prior, data_prior_file)
    geometry_error <- NULL
    tryCatch({
      marked_unc_names <- .pdmpsamplers_julia_call(
        "r_get_param_unc_names",
        normalizePath(stan_file, mustWork = TRUE),
        normalizePath(data_full_file, mustWork = TRUE)
      )
      marked_geometry <- build_marked_predictor_geometry(
        sdata, marked_unc_names, eligibility
      )
      marked_multipliers <- marked_observation_multipliers(
        sdata, eligibility$family
      )
    }, error = function(e) geometry_error <<- conditionMessage(e))
    if (!is.null(geometry_error)) {
      cli::cli_inform(c(
        "The affine marked geometry could not be certified; using the exact full-gradient sampler.",
        "i" = geometry_error
      ))
      subsampled <- FALSE
      use_hcv <- FALSE
      use_anchor_bank <- FALSE
      if (!flow %in% c("AdaptiveBoomerang", "PreconditionedZigZag",
                       "PreconditionedBPS", "DensePreconditionedZigZag",
                       "DensePreconditionedBPS")) {
        t_warmup <- requested_t_warmup
      }
      data_file <- tempfile(fileext = ".json")
      write_stan_json(sdata, data_file)
    } else if (identical(algorithm, "ThinningStrategy") &&
               !flow %in% c("Boomerang", "AdaptiveBoomerang") &&
               identical(eligibility$family, "gaussian") &&
               length(marked_geometry$designs) > 1L) {
      cli::cli_inform(c(
        "Subsampling is not eligible for this model; using the full-gradient sampler.",
        "i" = paste0(
          "ThinningStrategy has only an affine global clock, which cannot dominate ",
          "a distributional Gaussian residual envelope on a linear trajectory; ",
          "use GridThinningStrategy"
        )
      ))
      subsampled <- FALSE
      use_hcv <- FALSE
      use_anchor_bank <- FALSE
      if (!flow %in% c("PreconditionedZigZag", "PreconditionedBPS",
                       "DensePreconditionedZigZag", "DensePreconditionedBPS")) {
        t_warmup <- requested_t_warmup
      }
      data_file <- tempfile(fileext = ".json")
      write_stan_json(sdata, data_file)
    }
  } else {
    data_file <- tempfile(fileext = ".json")
    write_stan_json(sdata, data_file)
    if (!is.null(slab_prior)) {
      data_prior_file <- tempfile(fileext = ".json")
      write_stan_json(sdata_prior, data_prior_file)
    }
  }

  if (!is.null(save_model))
    cat(scode, file = save_model)

  csv_file <- tempfile(fileext = ".csv")
  jl_flow_mean <- if (is.null(flow_mean)) numeric(0) else flow_mean
  jl_flow_cov  <- if (is.null(flow_cov)) matrix(numeric(0), nrow = 0, ncol = 0) else flow_cov
  jl_discretize_dt <- if (is.null(discretize_dt)) 0.0 else discretize_dt

  # Validate sticky arguments (requires param_unc_names from BridgeStan)
  if (isTRUE(sticky)) {
    brms_prior <- as.data.frame(brms::prior_summary(empty_fit), stringsAsFactors = FALSE)
    fe_names <- unique(brms_prior$coef[brms_prior$class == "b" & nzchar(brms_prior$coef)])
    supported_coef_names <- supported_b_coef_names(
      fe_names = fe_names,
      formula = formula,
      data = data
    )

    data_for_names <- if (subsampled) data_full_file else data_file
    unc_names <- .pdmpsamplers_julia_call(
      "r_get_param_unc_names",
      normalizePath(stan_file, mustWork = TRUE),
      normalizePath(data_for_names, mustWork = TRUE)
    )
    if (!is.null(slab_prior)) {
      if (!is.null(kappa)) {
        cli::cli_abort("Use either legacy {.arg kappa} or dependent {.arg slab_prior}, not both.")
      }
      if (!is.slab_prior(slab_prior)) {
        cli::cli_abort("Argument {.arg slab_prior} must be created by a dependent slab constructor.")
      }
      if (is.null(model_prior) || !is.model_prior(model_prior)) {
        cli::cli_abort("Argument {.arg model_prior} must be provided when {.arg sticky} is {.code TRUE}.")
      }
      can_stick_full <- map_can_stick(
        unc_names,
        supported_coef_names = supported_coef_names,
        user_can_stick = can_stick
      )
      .validate_slab_prior_dimensions(slab_prior, length(unc_names), can_stick_full,
                                      unc_names = unc_names, model_prior = model_prior)
      sticky_args <- list(
        sticky = TRUE,
        can_stick = can_stick_full,
        model_prior = model_prior,
        parameter_prior = NULL,
        slab_prior = slab_prior,
        unc_names = unc_names,
        supported_coef_names = supported_coef_names
      )
    } else {
      sticky_args <- validate_brms_sticky(
        sticky, can_stick, model_prior, kappa,
        d = length(unc_names), unc_names = unc_names,
        supported_coef_names = supported_coef_names,
        prior = brms_prior, subsampled = subsampled
      )
      sticky_args$slab_prior <- NULL
      sticky_args$unc_names <- unc_names
    }
  } else {
    sticky_args <- list(sticky = FALSE, can_stick = NULL,
                        model_prior = NULL, parameter_prior = NULL,
                        slab_prior = NULL, unc_names = character(0))
  }

  if (subsampled) {
    jl_result <- .pdmpsamplers_julia_call(
      "r_pdmp_brms_marked",
      normalizePath(stan_file, mustWork = TRUE),
      normalizePath(data_full_file, mustWork = TRUE),
      normalizePath(data_prior_file, mustWork = TRUE),
      as.integer(N), subsample_size,
      eligibility$family, marked_geometry$designs,
      marked_geometry$design_indices, as.integer(marked_geometry$dimension),
      marked_geometry$offsets,
      marked_geometry$response, marked_geometry$se, marked_multipliers,
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
      seed = seed,
      compute_lp = compute_lp,
      use_fd_hvp = use_fd_hvp,
      post_warmup_simplify = post_warmup_simplify,
      sticky = sticky_args$sticky,
      can_stick = sticky_args$can_stick,
      model_prior = sticky_args$model_prior,
      parameter_prior = sticky_args$parameter_prior,
      n_anchor_updates = as.integer(n_anchor_updates),
      use_anchor_bank = isTRUE(use_anchor_bank),
      bank_capacity = as.integer(bank_capacity),
      use_hcv = isTRUE(use_hcv)
    )
  } else {
    jl_result <- .pdmpsamplers_julia_call(
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
      seed = seed,
      compute_lp = compute_lp,
      use_fd_hvp = use_fd_hvp,
      post_warmup_simplify = post_warmup_simplify,
      sticky = sticky_args$sticky,
      can_stick = sticky_args$can_stick,
      model_prior = sticky_args$model_prior,
      parameter_prior = sticky_args$parameter_prior,
      slab_prior = sticky_args$slab_prior,
      prior_data_file = if (!is.null(sticky_args$slab_prior)) normalizePath(data_prior_file, mustWork = TRUE) else NULL,
      unc_names = sticky_args$unc_names
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
  attr(empty_fit, "bridge_call_counts") <- jl_result$bridge_call_counts
  attr(empty_fit, "marked_subsampling") <- isTRUE(jl_result$marked_subsampling)
  class(empty_fit) <- unique(c("pdmp_brmsfit", class(empty_fit)))

  if (isTRUE(sticky_args$sticky) && !is.null(jl_result$inclusion_probs)) {
    incl_list <- build_sticky_inclusion_probs(
      incl_raw = jl_result$inclusion_probs,
      can_stick = sticky_args$can_stick,
      unc_names = unc_names,
      supported_coef_names = sticky_args$supported_coef_names
    )

    attr(empty_fit, "sticky") <- list(
      inclusion_probs = incl_list,
      can_stick = sticky_args$can_stick,
      supported_coef_names = sticky_args$supported_coef_names,
      unc_names = unc_names,
      model_prior = sticky_args$model_prior
    )
    class(empty_fit) <- unique(c("sticky_brmsfit", class(empty_fit)))
  }

  empty_fit
}

build_sticky_inclusion_probs <- function(incl_raw, can_stick, unc_names, supported_coef_names = NULL) {
  if (is.environment(incl_raw))
    incl_raw <- as.list(incl_raw)
  if (is.numeric(incl_raw))
    incl_raw <- list(chain1 = incl_raw)
  if (!is.list(incl_raw))
    cli::cli_abort("Sticky inclusion probabilities must be a numeric vector or a list of numeric vectors.")

  stickable_idx <- which(can_stick)
  if (!is.null(supported_coef_names)) {
    resolved <- .resolve_stickable_mapping(
      unc_names,
      supported_coef_names = supported_coef_names
    )
    name_lookup <- stats::setNames(resolved$display_names, as.character(resolved$indices))
    stickable_names <- unname(name_lookup[as.character(stickable_idx)])
    missing_names <- is.na(stickable_names)
    if (any(missing_names))
      stickable_names[missing_names] <- unc_names[stickable_idx[missing_names]]
  } else {
    stickable_names <- unc_names[stickable_idx]
  }

  incl_list <- lapply(incl_raw, function(probs) {
    if (!is.numeric(probs))
      cli::cli_abort("Each chain's sticky inclusion probabilities must be numeric.")
    if (length(probs) != length(unc_names))
      cli::cli_abort(c(
        "Sticky inclusion probabilities had unexpected length.",
        "x" = "Expected length {length(unc_names)} and got {length(probs)}."
      ))

    named_probs <- probs[stickable_idx]
    named_probs[!is.finite(named_probs)] <- 0.0
    named_probs <- pmin(pmax(named_probs, 0.0), 1.0)
    names(named_probs) <- stickable_names
    named_probs
  })

  names(incl_list) <- paste0("chain", seq_along(incl_list))

  incl_list
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
