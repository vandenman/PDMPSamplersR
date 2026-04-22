# Integration tests for brm_pdmp() variable selection (sticky sampling).
#
# These tests run full PDMP sampling with sticky dynamics and verify
# inclusion probabilities recover ground truth.
# Gated by PDMPSAMPLERSR_SLOW_TESTS + Julia/BridgeStan availability.
#
# Run locally:
#   NOT_CRAN=true PDMPSAMPLERSR_SLOW_TESTS=true Rscript -e 'require(PDMPSamplersR); testthat::test_file("tests/testthat/test-brms_varsel.R")'

skip_if_no_slow_tests <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not(
    identical(Sys.getenv("PDMPSAMPLERSR_SLOW_TESTS"), "true"),
    "Slow MCMC tests not enabled (set PDMPSAMPLERSR_SLOW_TESTS=true)"
  )
  skip_if_no_brms_setup()
}

# ==============================================================================
# Shared fixture infrastructure
# ==============================================================================
#
# Expensive brm_pdmp() fits are cached and reused across test blocks to
# keep total file runtime under ~300 s.

.varsel_cache <- new.env(parent = emptyenv())

shared_gaussian_data <- function() {
  if (!is.null(.varsel_cache$df)) return(.varsel_cache$df)
  set.seed(42)
  n <- 100L; d <- 6L
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(3.0, -2.0, rep(0, d - 2))
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  .varsel_cache$d <- d
  .varsel_cache$df <- data.frame(y = as.numeric(y), X)
  .varsel_cache$df
}

shared_gaussian_sticky_fit <- function() {
  if (!is.null(.varsel_cache$fit)) return(.varsel_cache$fit)
  df <- shared_gaussian_data()
  .varsel_cache$fit <- brm_pdmp(
    y ~ x1 + x2 + x3 + x4 + x5 + x6, data = df, family = gaussian(),
    prior = brms::prior(normal(0, 5), class = b),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 10000, show_progress = FALSE,
    sticky = TRUE, model_prior = bernoulli(0.5)
  )
  .varsel_cache$fit
}

# ==============================================================================
# Section 1: Sparse Gaussian regression
# ==============================================================================

test_that("brm_pdmp sticky: sparse Gaussian regression recovers active predictors", {
  skip_if_no_slow_tests()

  fit <- shared_gaussian_sticky_fit()
  d <- .varsel_cache$d

  expect_s3_class(fit, "brmsfit")

  sticky_meta <- attr(fit, "sticky")
  expect_false(is.null(sticky_meta))
  expect_false(is.null(sticky_meta$inclusion_probs))

  incl <- sticky_meta$inclusion_probs$chain1
  expect_length(incl, d)

  active_incl <- incl[1:2]
  inactive_incl <- incl[3:d]
  expect_true(all(active_incl > 0.5),
              info = paste("Active inclusion probs:", paste(round(active_incl, 3), collapse = ", ")))
  expect_true(all(inactive_incl < 0.5),
              info = paste("Inactive inclusion probs:", paste(round(inactive_incl, 3), collapse = ", ")))

  expect_true(min(active_incl) > max(inactive_incl))
})

# ==============================================================================
# Section 2: Sparse logistic regression
# ==============================================================================

test_that("brm_pdmp sticky: sparse logistic regression", {
  skip_if_no_slow_tests()

  set.seed(123)
  n <- 200
  d <- 6
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(2.0, -1.5, rep(0, d - 2))
  prob <- plogis(X %*% beta_true)
  y <- rbinom(n, 1, prob)
  df <- data.frame(y = y, X)

  formula <- y ~ x1 + x2 + x3 + x4 + x5 + x6
  fit <- brm_pdmp(
    formula, data = df, family = brms::bernoulli(),
    prior = brms::prior(normal(0, 5), class = b),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 10000, show_progress = FALSE,
    sticky = TRUE, model_prior = bernoulli(0.5)
  )

  expect_s3_class(fit, "brmsfit")
  sticky_meta <- attr(fit, "sticky")
  expect_false(is.null(sticky_meta))

  incl <- sticky_meta$inclusion_probs$chain1
  active_incl <- incl[1:2]
  inactive_incl <- incl[3:d]
  expect_true(all(active_incl > 0.5),
              info = paste("Active:", paste(round(active_incl, 3), collapse = ", ")))
  expect_true(min(active_incl) > max(inactive_incl))
})

# ==============================================================================
# Section 3: Intercept handling + default non-sticky
# ==============================================================================

test_that("brm_pdmp sticky: intercept is never stickable by default", {
  skip_if_no_slow_tests()

  fit <- shared_gaussian_sticky_fit()
  sticky_meta <- attr(fit, "sticky")
  expect_false(is.null(sticky_meta))

  intercept_idx <- grep("Intercept", sticky_meta$unc_names)
  expect_true(length(intercept_idx) > 0)
  expect_false(any(sticky_meta$can_stick[intercept_idx]))
})

test_that("brm_pdmp: sticky=FALSE (default) unchanged", {
  skip_if_no_slow_tests()

  df <- shared_gaussian_data()
  fit <- brm_pdmp(y ~ x1 + x2 + x3 + x4 + x5 + x6, data = df,
                  family = gaussian(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 3000, show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")
  expect_null(attr(fit, "sticky"))
})

# ==============================================================================
# Section 4: Error cases
# ==============================================================================

test_that("brm_pdmp sticky: errors without model_prior", {
  skip_if_no_slow_tests()

  df <- data.frame(y = rnorm(10), x = rnorm(10))
  expect_error(
    brm_pdmp(y ~ x, data = df, family = gaussian(),
             prior = brms::prior(normal(0, 5), class = b),
             flow = "ZigZag", algorithm = "GridThinningStrategy",
             T = 1000, show_progress = FALSE, sticky = TRUE),
    "model_prior"
  )
})

test_that("brm_pdmp sticky: unsupported prior triggers informative error", {
  skip_if_no_slow_tests()

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  expect_error(
    brm_pdmp(y ~ x, data = df, family = gaussian(),
             prior = brms::prior(cauchy(0, 2.5), class = b),
             flow = "ZigZag", algorithm = "GridThinningStrategy",
             T = 1000, show_progress = FALSE,
             sticky = TRUE, model_prior = bernoulli(0.5)),
    "allowlist|parameter_prior"
  )
})

# ==============================================================================
# Section 5: Subsampled sticky gradients
# ==============================================================================

test_that("brm_pdmp sticky: works with subsampled gradients", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 200
  d <- 6
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(3.0, -2.0, rep(0, d - 2))
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = as.numeric(y), X)

  formula <- y ~ x1 + x2 + x3 + x4 + x5 + x6
  fit <- brm_pdmp(
    formula, data = df, family = gaussian(),
    prior = brms::set_prior("normal(0, 5)", class = "b"),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 10000, t_warmup = 2000, show_progress = FALSE,
    sticky = TRUE, model_prior = bernoulli(0.5),
    subsample_size = 50L
  )

  expect_s3_class(fit, "brmsfit")
  sticky_meta <- attr(fit, "sticky")
  expect_false(is.null(sticky_meta))
  expect_false(is.null(sticky_meta$inclusion_probs))

  incl <- sticky_meta$inclusion_probs$chain1
  expect_length(incl, d)

  active_incl <- incl[1:2]
  inactive_incl <- incl[3:d]
  expect_true(min(active_incl) > max(inactive_incl),
              info = paste("Active:", paste(round(active_incl, 3), collapse = ", "),
                           "; max inactive:", round(max(inactive_incl), 3)))
})

# ==============================================================================
# Section 6: Beta-Bernoulli model prior
# ==============================================================================

test_that("brm_pdmp sticky: beta-bernoulli prior works", {
  skip_if_no_slow_tests()

  df <- shared_gaussian_data()
  d <- .varsel_cache$d

  fit <- brm_pdmp(
    y ~ x1 + x2 + x3 + x4 + x5 + x6, data = df, family = gaussian(),
    prior = brms::prior(normal(0, 5), class = b),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 8000, show_progress = FALSE,
    sticky = TRUE, model_prior = betabernoulli(1, 1)
  )

  sticky_meta <- attr(fit, "sticky")
  expect_false(is.null(sticky_meta))

  incl <- sticky_meta$inclusion_probs$chain1
  active_incl <- incl[1:2]
  inactive_incl <- incl[3:d]
  expect_true(min(active_incl) > max(inactive_incl))
})

# ==============================================================================
# Section 7: Correlated design
# ==============================================================================

test_that("brm_pdmp sticky: correlated design surfaces inclusion dilution", {
  skip_if_no_slow_tests()

  set.seed(314)
  n <- 200
  d <- 6
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  x1 <- z1
  x2 <- 0.8 * z1 + 0.6 * z2  # cor(x1, x2) ~ 0.8
  X <- cbind(x1, x2, matrix(rnorm(n * (d - 2)), n, d - 2))
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(2.0, 0, 1.5, 0, 0, 0)
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = as.numeric(y), X)

  formula <- y ~ x1 + x2 + x3 + x4 + x5 + x6
  fit <- brm_pdmp(
    formula, data = df, family = gaussian(),
    prior = brms::set_prior("normal(0, 5)", class = "b"),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 15000, show_progress = FALSE,
    sticky = TRUE, model_prior = bernoulli(0.5)
  )

  sticky_meta <- attr(fit, "sticky")
  expect_false(is.null(sticky_meta))
  incl <- sticky_meta$inclusion_probs$chain1
  expect_length(incl, d)

  expect_gt(incl[3], 0.7,
            label = paste("x3 inclusion:", round(incl[3], 3)))

  expect_lt(incl[4], 0.3,
            label = paste("x4 inclusion:", round(incl[4], 3)))
  expect_lt(incl[5], 0.3,
            label = paste("x5 inclusion:", round(incl[5], 3)))

  expect_gt(incl[1] + incl[2], 0.8,
            label = paste("x1+x2 combined:", round(incl[1] + incl[2], 3)))

  draws <- posterior::as_draws_matrix(fit)
  b_cols <- sort(grep("^b_(?!Intercept)", colnames(draws), perl = TRUE, value = TRUE))
  b_draws <- as.matrix(draws[, b_cols])
  model_indicators <- apply(1L * (b_draws != 0), 1, paste0, collapse = "")
  model_freq <- sort(table(model_indicators), decreasing = TRUE)
  model_probs <- model_freq / sum(model_freq)

  true_model <- "101000"
  competitor <- "011000"

  top_models <- names(model_probs)[1:min(5, length(model_probs))]
  expect_true(
    true_model %in% top_models || competitor %in% top_models,
    info = paste("Top 5 models:", paste(top_models, collapse = ", "),
                 "do not contain true or competitor model")
  )

  p_true <- if (true_model %in% names(model_probs)) model_probs[true_model] else 0
  p_comp <- if (competitor %in% names(model_probs)) model_probs[competitor] else 0
  expect_gte(as.numeric(p_true), as.numeric(p_comp) * 0.5,
             label = paste("P(true):", round(p_true, 4),
                           "P(competitor):", round(p_comp, 4)))
})

# ==============================================================================
# Section 8: Anchor-update invariance under subsampled sticky
# ==============================================================================

test_that("brm_pdmp sticky subsampled: anchor-update count does not alter active/inactive ranking", {
  skip_if_no_slow_tests()

  set.seed(7)
  n <- 200
  d <- 6
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(3.0, -2.0, rep(0, d - 2))
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = as.numeric(y), X)
  formula <- y ~ x1 + x2 + x3 + x4 + x5 + x6

  fit_few <- brm_pdmp(
    formula, data = df, family = gaussian(),
    prior = brms::set_prior("normal(0, 5)", class = "b"),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 10000, t_warmup = 2000, show_progress = FALSE,
    sticky = TRUE, model_prior = bernoulli(0.5),
    subsample_size = 50L, n_anchor_updates = 3L
  )
  fit_many <- brm_pdmp(
    formula, data = df, family = gaussian(),
    prior = brms::set_prior("normal(0, 5)", class = "b"),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 10000, t_warmup = 2000, show_progress = FALSE,
    sticky = TRUE, model_prior = bernoulli(0.5),
    subsample_size = 50L, n_anchor_updates = 15L
  )

  incl_few  <- attr(fit_few,  "sticky")$inclusion_probs$chain1
  incl_many <- attr(fit_many, "sticky")$inclusion_probs$chain1

  expect_true(min(incl_few[1:2])  > max(incl_few[3:d]),
              info = paste("n_anchor=3 ; active:", paste(round(incl_few[1:2], 3), collapse = ", "),
                           "; max inactive:", round(max(incl_few[3:d]), 3)))
  expect_true(min(incl_many[1:2]) > max(incl_many[3:d]),
              info = paste("n_anchor=15; active:", paste(round(incl_many[1:2], 3), collapse = ", "),
                           "; max inactive:", round(max(incl_many[3:d]), 3)))
})

# ==============================================================================
# Section 9: Summary helpers (inclusion_prob, median_probability_model,
#            model_averaged_mean)
# ==============================================================================

test_that("inclusion_prob returns named vector in [0,1] on a real sticky fit", {
  skip_if_no_slow_tests()

  fit <- shared_gaussian_sticky_fit()
  d <- .varsel_cache$d

  ip <- inclusion_prob(fit)
  expect_true(is.numeric(ip))
  expect_length(ip, d)
  expect_true(all(ip >= 0 & ip <= 1))
  expect_named(ip)
})

test_that("median_probability_model recovers active predictors on a real sticky fit", {
  skip_if_no_slow_tests()

  fit <- shared_gaussian_sticky_fit()

  expect_warning(
    mpm <- median_probability_model(fit),
    "collinearity"
  )
  expect_true(is.character(mpm))
  expect_true("b.x1" %in% mpm)
  expect_true("b.x2" %in% mpm)
  expect_false("b.x3" %in% mpm)
})

test_that("model_averaged_mean returns named numeric for real sticky fit", {
  skip_if_no_slow_tests()

  fit <- shared_gaussian_sticky_fit()

  mam <- model_averaged_mean(fit)
  expect_true(is.numeric(mam))
  expect_true(length(mam) > 0L)
  expect_named(mam)
  expect_true(abs(mam["x1"]) > abs(mam["x3"]),
              info = "Model-averaged mean for active predictor should dominate inactive one")
})
