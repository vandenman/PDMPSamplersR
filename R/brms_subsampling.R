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
  Y_sub <- sdata$Y[indices]
  if (length(Y_sub) == 1L) Y_sub <- array(Y_sub)
  sub$Y <- Y_sub
  sub$X <- sdata$X[indices, , drop = FALSE]
  sub$means_X <- means_X
  sub
}

make_prior_standata <- function(sdata, means_X) {
  K <- sdata$K
  prior <- list(
    N = 1L,
    Y = array(if (is.integer(sdata$Y)) 0L else 0.0),
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
#' @description
#' Deprecated: use [brm_pdmp()] with the `subsample_size` argument instead.
#'
#' @inheritParams brm_pdmp
#' @param subsample_size Integer, number of observations per subsample.
#'   Default: `max(1, nrow(data) %/% 10)`.
#' @param n_anchor_updates Integer, number of anchor updates during warmup
#'   (default: 10).
#'
#' @return A `brmsfit` object.
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
  cli::cli_warn(c(
    "{.fn brm_pdmp_subsampled} is deprecated.",
    "i" = "Use {.fn brm_pdmp} with {.arg subsample_size} instead."
  ))

  if (is.null(subsample_size))
    subsample_size <- max(1L, nrow(data) %/% 10L)

  brm_pdmp(
    formula = formula, data = data, family = family,
    prior = prior,
    flow = flow, algorithm = algorithm,
    adaptive_scheme = adaptive_scheme,
    T = T, t0 = t0, t_warmup = t_warmup,
    flow_mean = flow_mean, flow_cov = flow_cov, c0 = c0,
    grid_n = grid_n, grid_t_max = grid_t_max,
    show_progress = show_progress,
    discretize_dt = discretize_dt,
    n_chains = n_chains, threaded = threaded,
    compute_lp = compute_lp,
    subsample_size = subsample_size,
    n_anchor_updates = n_anchor_updates,
    stanvars = stanvars, sample_prior = sample_prior,
    save_model = save_model,
    ...
  )
}
