fix_brms_stancode <- function(stancode) {
  lines <- strsplit(stancode, "\n")[[1]]
  n <- length(lines)

  td_start <- grep("^transformed data\\s*\\{", lines)
  if (length(td_start) == 0) return(stancode)

  means_decl <- grep("^\\s*vector\\[Kc\\] means_X;", lines)
  means_assign <- grep("means_X\\[i - 1\\] = mean\\(X\\[, i\\]\\);", lines)
  if (length(means_decl) == 0 || length(means_assign) == 0) return(stancode)

  data_close <- grep("^\\s*int prior_only;", lines)
  if (length(data_close) == 0) {
    data_close <- grep("^\\}", lines)
    data_close <- data_close[data_close < td_start[1]]
    if (length(data_close) == 0) return(stancode)
    data_close <- max(data_close)
  } else {
    data_close <- data_close[1]
  }

  new_lines <- character(0)
  skip <- integer(0)
  skip <- c(means_decl[1], means_assign[1])
  inserted <- FALSE

  for (i in seq_len(n)) {
    if (i %in% skip) next
    new_lines <- c(new_lines, lines[i])
    if (i == data_close && !inserted) {
      new_lines <- c(new_lines, "  vector[Kc] means_X;  // column means of X (fixed for subsampling)")
      inserted <- TRUE
    }
  }

  paste(new_lines, collapse = "\n")
}

subset_standata <- function(sdata, indices, means_X) {
  sub <- sdata
  sub$N <- length(indices)
  sub$Y <- sdata$Y[indices]
  sub$X <- sdata$X[indices, , drop = FALSE]
  sub$means_X <- means_X
  sub
}

make_prior_standata <- function(sdata, means_X) {
  K <- sdata$K
  prior <- list(
    N = 1L,
    Y = if (is.integer(sdata$Y)) 0L else 0.0,
    K = K,
    Kc = sdata$Kc,
    X = matrix(0, nrow = 1, ncol = K),
    means_X = means_X,
    prior_only = 1L
  )
  prior
}

#' Fit a brms model using PDMP samplers with subsampled gradients
#'
#' Uses BridgeStan data-swapping to compute control-variate subsampled
#' gradients. The same compiled Stan model is shared across the full dataset,
#' subsample, and prior-only evaluations.
#'
#' @inheritParams brm_pdmp
#' @param subsample_size Integer, number of observations per subsample.
#'   Must be less than `nrow(data)`. Default: `max(1, nrow(data) %/% 10)`.
#' @param n_anchor_updates Integer, number of anchor updates during warmup
#'   (default: 10).
#'
#' @return A `brmsfit` object.
#'
#' @details
#' The subsampled gradient uses a control-variate estimator that combines
#' a subsample gradient, a prior-only gradient, and a cached full-data
#' correction term (updated at anchor points during warmup).
#'
#' A centering fix is applied to the brms-generated Stan code to ensure
#' consistent predictor centering across subsamples. This moves `means_X`
#' from `transformed data` to `data`, always using the full-data column
#' means.
#'
#' Currently supports fixed-effects models only (no random effects).
#'
#' @export
brm_pdmp_subsampled <- function(
    formula, data, family = gaussian(),
    prior = NULL,
    subsample_size = NULL,
    flow = c("ZigZag", "BouncyParticle", "Boomerang",
             "AdaptiveBoomerang", "PreconditionedZigZag", "PreconditionedBPS"),
    algorithm = c("GridThinningStrategy", "ThinningStrategy",
                  "RootsPoissonStrategy"),
    adaptive_scheme = c("diagonal", "fullrank"),
    T = 50000, t0 = 0.0, t_warmup = 1000.0,
    n_anchor_updates = 10L,
    flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
    grid_n = 30, grid_t_max = 2.0,
    show_progress = TRUE,
    discretize_dt = NULL,
    n_chains = 1L, threaded = FALSE,
    compute_lp = FALSE,
    stanvars = NULL, sample_prior = "no",
    save_model = NULL,
    ...
) {
  if (!requireNamespace("brms", quietly = TRUE))
    cli::cli_abort("Package {.pkg brms} is required for {.fn brm_pdmp_subsampled}.")
  if (!requireNamespace("rstan", quietly = TRUE))
    cli::cli_abort("Package {.pkg rstan} is required for {.fn brm_pdmp_subsampled}.")

  if (sample_prior != "no")
    cli::cli_abort('{.fn brm_pdmp_subsampled} only supports {.code sample_prior = "no"}.')

  flow <- match.arg(flow)
  algorithm <- match.arg(algorithm)
  adaptive_scheme <- match.arg(adaptive_scheme)
  N <- nrow(data)

  if (is.null(subsample_size)) subsample_size <- max(1L, N %/% 10L)
  subsample_size <- as.integer(subsample_size)
  if (subsample_size >= N)
    cli::cli_abort("{.arg subsample_size} ({subsample_size}) must be less than nrow(data) ({N}).")

  n_anchor_updates <- as.integer(n_anchor_updates)

  if (flow == "AdaptiveBoomerang") {
    if (algorithm != "GridThinningStrategy")
      cli::cli_abort("{.val AdaptiveBoomerang} requires {.val GridThinningStrategy}.")
    if (t_warmup == 0) {
      t_warmup <- (T - t0) / 5
      cli::cli_inform("Setting {.arg t_warmup} to {t_warmup} (20% of sampling time) for {.val AdaptiveBoomerang}.")
    }
  }
  if (flow %in% c("PreconditionedZigZag", "PreconditionedBPS")) {
    if (algorithm != "GridThinningStrategy")
      cli::cli_abort("{.val {flow}} requires {.val GridThinningStrategy}.")
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

  scode_fixed <- fix_brms_stancode(scode)

  means_X <- colMeans(sdata$X[, -1, drop = FALSE])

  sdata_full <- sdata
  sdata_full$means_X <- means_X

  sdata_prior <- make_prior_standata(sdata, means_X)

  init_indices <- sample.int(N, subsample_size)
  sdata_sub <- subset_standata(sdata, init_indices, means_X)

  empty_fit <- brms::brm(formula, data = data, family = family,
                         prior = prior, stanvars = stanvars,
                         sample_prior = sample_prior,
                         empty = TRUE, ...)

  stan_file <- tempfile(fileext = ".stan")
  cat(scode_fixed, file = stan_file)
  data_full_file <- tempfile(fileext = ".json")
  write_stan_json(sdata_full, data_full_file)
  data_prior_file <- tempfile(fileext = ".json")
  write_stan_json(sdata_prior, data_prior_file)
  data_sub_file <- tempfile(fileext = ".json")
  write_stan_json(sdata_sub, data_sub_file)

  if (!is.null(save_model))
    cat(scode_fixed, file = save_model)

  csv_file <- tempfile(fileext = ".csv")

  jl_flow_mean <- if (is.null(flow_mean)) numeric(0) else flow_mean
  jl_flow_cov  <- if (is.null(flow_cov)) matrix(numeric(0), nrow = 0, ncol = 0) else flow_cov
  jl_discretize_dt <- if (is.null(discretize_dt)) 0.0 else discretize_dt

  # Serialize full standata fields needed for subsample JSON construction
  Y_full <- as.numeric(sdata$Y)
  X_full <- sdata$X

  csv_paths <- JuliaCall::julia_call(
    "r_pdmp_brms_subsampled",
    normalizePath(stan_file, mustWork = TRUE),
    normalizePath(data_full_file, mustWork = TRUE),
    normalizePath(data_prior_file, mustWork = TRUE),
    normalizePath(data_sub_file, mustWork = TRUE),
    Y_full, X_full, means_X,
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
    compute_lp = compute_lp
  )

  stanfit <- rstan::read_stan_csv(csv_paths)
  empty_fit$fit <- stanfit
  empty_fit <- brms::rename_pars(empty_fit)
  empty_fit
}
