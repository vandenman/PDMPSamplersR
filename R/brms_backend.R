random_effect_subsampling_diagnostics <- function(sdata, subsample_size) {
  N <- sdata$N %||% 0L
  if (!is.numeric(subsample_size) || length(subsample_size) != 1L ||
      !is.numeric(N) || length(N) != 1L || subsample_size <= 0L || subsample_size >= N) {
    return(data.frame())
  }

  z_names <- grep("^Z_[0-9]+_[0-9]+$", names(sdata), value = TRUE)
  if (length(z_names) == 0L) return(data.frame())

  out <- lapply(z_names, function(name) {
    Z <- sdata[[name]]
    if (!is.matrix(Z) || nrow(Z) != N || ncol(Z) == 0L) return(NULL)

    support <- colSums(abs(Z) > 0)
    support <- support[is.finite(support) & support > 0]
    if (length(support) == 0L) return(NULL)

    min_support <- min(support)
    data.frame(
      block = name,
      n_columns = ncol(Z),
      min_support = min_support,
      expected_support = subsample_size * min_support / N,
      p_zero = stats::dhyper(0L, min_support, N - min_support, subsample_size)
    )
  })

  out <- Filter(Negate(is.null), out)
  if (length(out) == 0L) return(data.frame())
  do.call(rbind, out)
}

warn_if_low_random_effect_subsampling_support <- function(sdata, subsample_size) {
  diag <- random_effect_subsampling_diagnostics(sdata, subsample_size)
  if (!nrow(diag)) return(invisible(NULL))

  worst_idx <- which.min(diag$expected_support)
  worst <- diag[worst_idx, , drop = FALSE]
  if (worst$expected_support[[1]] >= 5 && worst$p_zero[[1]] <= 0.01) {
    return(invisible(NULL))
  }

  cli::cli_warn(c(
    "Subsampled random-effects gradients may be unstable for this model.",
    "x" = "The sparsest random-effects design block {.code {worst$block[[1]]}} has minimum support {worst$min_support[[1]]} observations per coefficient in the full data, but only an expected {formatC(worst$expected_support[[1]], digits = 2, format = 'f')} observations per minibatch at {.arg subsample_size} = {subsample_size}.",
    "i" = "A coefficient in that block is absent from a minibatch with probability about {formatC(worst$p_zero[[1]], digits = 3, format = 'f')}.",
    "i" = "This mainly affects group-level parameters: fixed effects can look stable while random effects still drift.",
    "i" = "Consider increasing {.arg subsample_size} or using full-data gradients for hierarchical models."
  ))

  invisible(NULL)
}

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
#' @param seed NULL (default) or a non-negative integer seed passed through to
#'   Julia's sampler RNG.
#' @param compute_lp Logical; compute `lp__` via `BridgeStan::log_density()`
#'   for each sample (default: FALSE). Adds overhead but enables
#'   `bridge_sampler()` and populates the `lp__` diagnostic column.
#' @param subsample_size Integer number of observations per subsample,
#'   or NULL (default) for full-data gradients. When non-NULL, a
#'   BridgeStan control-variate subsampled gradient is used. Must be
#'   less than `nrow(data)`. Fixed-effects and random-effects models are
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
#' @param sticky Logical; enable spike-and-slab variable selection for
#'   population-level coefficients (default: FALSE). Requires
#'   `model_prior` to be set.
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
    n_chains = 1L, threaded = FALSE, seed = NULL,
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
  hvp_mode <- match.arg(hvp_mode)

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
  if (subsampled && !is.null(slab_prior)) {
    cli::cli_abort("Dependent {.arg slab_prior} is not yet supported together with {.arg subsample_size}.")
  }
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
    warn_if_low_random_effect_subsampling_support(sdata, subsample_size)
  }

  if (subsampled) {
    scode_ext <- make_ext_cpp_stancode(scode, formula, data, family,
                                       prior, stanvars, sample_prior, ...)
    Y_full <- as.numeric(sdata$Y)
    X_full <- sdata$X
  }
  if (subsampled || !is.null(slab_prior)) {
    sdata_prior <- make_prior_standata(sdata)
  }

  empty_fit <- brms::brm(formula, data = data, family = family,
                         prior = prior, stanvars = stanvars,
                         sample_prior = sample_prior,
                         empty = TRUE, ...)

  stan_file <- cached_stan_model(scode)

  if (subsampled) {
    stan_file_ext <- cached_stan_model(scode_ext)
    data_full_file <- tempfile(fileext = ".json")
    write_stan_json(sdata, data_full_file)
    data_prior_file <- tempfile(fileext = ".json")
    write_stan_json(sdata_prior, data_prior_file)
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
  jl_resample_dt <- if (is.null(resample_dt)) 0.0 else resample_dt

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
      seed = seed,
      compute_lp = compute_lp,
      resample_dt = jl_resample_dt,
      hvp_mode = hvp_mode,
      use_hcv = use_hcv,
      use_anchor_bank = use_anchor_bank,
      bank_capacity = as.integer(bank_capacity),
      use_fd_hvp = use_fd_hvp,
      post_warmup_simplify = post_warmup_simplify,
      use_fd_hcv = use_fd_hcv,
      sticky = sticky_args$sticky,
      can_stick = sticky_args$can_stick,
      model_prior = sticky_args$model_prior,
      parameter_prior = sticky_args$parameter_prior,
      slab_prior = sticky_args$slab_prior,
      unc_names = sticky_args$unc_names
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
