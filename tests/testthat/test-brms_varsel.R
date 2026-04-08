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
# Section 1: Sparse Gaussian regression
# ==============================================================================

test_that("brm_pdmp sticky: sparse Gaussian regression recovers active predictors", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 100
  d <- 10
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(2.0, -1.5, 1.0, rep(0, d - 3))
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = as.numeric(y), X)

  formula <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
  fit <- brm_pdmp(
    formula, data = df, family = gaussian(),
    prior = brms::prior(normal(0, 5), class = b),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 100000, show_progress = FALSE,
    sticky = TRUE, model_prior = bernoulli(0.5)
  )

  expect_s3_class(fit, "brmsfit")

  sticky_meta <- attr(fit, "sticky")
  expect_false(is.null(sticky_meta))
  expect_false(is.null(sticky_meta$inclusion_probs))

  incl <- sticky_meta$inclusion_probs$chain1
  expect_length(incl, d)

  active_incl <- incl[1:3]
  inactive_incl <- incl[4:d]
  expect_true(all(active_incl > 0.5),
              info = paste("Active inclusion probs:", paste(round(active_incl, 3), collapse = ", ")))
  expect_true(all(inactive_incl < 0.5),
              info = paste("Inactive inclusion probs:", paste(round(inactive_incl, 3), collapse = ", ")))

  # Active predictors should rank above inactive
  expect_true(min(active_incl) > max(inactive_incl))
})

# ==============================================================================
# Section 2: Sparse logistic regression
# ==============================================================================

test_that("brm_pdmp sticky: sparse logistic regression", {
  skip_if_no_slow_tests()

  set.seed(123)
  n <- 200
  d <- 8
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(1.5, -1.0, rep(0, d - 2))
  prob <- plogis(X %*% beta_true)
  y <- rbinom(n, 1, prob)
  df <- data.frame(y = y, X)

  formula <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
  fit <- brm_pdmp(
    formula, data = df, family = brms::bernoulli(),
    prior = brms::prior(normal(0, 5), class = b),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 100000, show_progress = FALSE,
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
# Section 3: Intercept handling
# ==============================================================================

test_that("brm_pdmp sticky: intercept is never stickable by default", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 2 + 0.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(
    y ~ x, data = df, family = gaussian(),
    prior = brms::prior(normal(0, 5), class = b),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 50000, show_progress = FALSE,
    sticky = TRUE, model_prior = bernoulli(0.5)
  )

  sticky_meta <- attr(fit, "sticky")
  expect_false(is.null(sticky_meta))

  # Intercept should not be stickable
  intercept_idx <- grep("Intercept", sticky_meta$unc_names)
  expect_true(length(intercept_idx) > 0)
  expect_false(any(sticky_meta$can_stick[intercept_idx]))
})

# ==============================================================================
# Section 4: Default non-sticky behavior unchanged
# ==============================================================================

test_that("brm_pdmp: sticky=FALSE (default) unchanged", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 2 + 0.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = gaussian(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")
  expect_null(attr(fit, "sticky"))
})

# ==============================================================================
# Section 5: Error cases
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

test_that("brm_pdmp sticky: errors with subsampling", {
  skip_if_no_slow_tests()

  df <- data.frame(y = rnorm(100), x = rnorm(100))
  expect_error(
    brm_pdmp(y ~ x, data = df, family = gaussian(),
             prior = brms::prior(normal(0, 5), class = b),
             flow = "ZigZag", algorithm = "GridThinningStrategy",
             T = 1000, show_progress = FALSE,
             sticky = TRUE, model_prior = bernoulli(0.5),
             subsample_size = 50L),
    "subsampled"
  )
})

# ==============================================================================
# Section 6: Beta-Bernoulli model prior
# ==============================================================================

test_that("brm_pdmp sticky: beta-bernoulli prior works", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 100
  d <- 6
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(2.0, -1.5, rep(0, d - 2))
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = as.numeric(y), X)

  formula <- y ~ x1 + x2 + x3 + x4 + x5 + x6
  fit <- brm_pdmp(
    formula, data = df, family = gaussian(),
    prior = brms::prior(normal(0, 5), class = b),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 100000, show_progress = FALSE,
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
# Section 7: Unsupported prior errors
# ==============================================================================

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
# Section 8: Correlated design
# ==============================================================================

test_that("brm_pdmp sticky: correlated design surfaces inclusion dilution", {
  skip_if_no_slow_tests()

  # Design: x1 and x2 are correlated (rho ~ 0.8), only x1 is truly active.
  # x3 is active and uncorrelated. x4, x5 are inactive.
  # Under correlation, marginal inclusion probability for x1 may be diluted
  # toward x2, while model probabilities reveal the correct grouping.
  set.seed(314)
  n <- 200
  d <- 5
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  x1 <- z1
  x2 <- 0.8 * z1 + 0.6 * z2  # cor(x1, x2) ~ 0.8
  X <- cbind(x1, x2, matrix(rnorm(n * (d - 2)), n, d - 2))
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(1.5, 0, 1.0, 0, 0)
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = as.numeric(y), X)

  formula <- y ~ x1 + x2 + x3 + x4 + x5
  fit <- brm_pdmp(
    formula, data = df, family = gaussian(),
    prior = brms::set_prior("normal(0, 5)", class = "b"),
    flow = "ZigZag", algorithm = "GridThinningStrategy",
    T = 200000, show_progress = FALSE,
    sticky = TRUE, model_prior = bernoulli(0.5)
  )

  sticky_meta <- attr(fit, "sticky")
  expect_false(is.null(sticky_meta))
  incl <- sticky_meta$inclusion_probs$chain1
  expect_length(incl, d)

  # Coefficients are in order: x1, x2, x3, x4, x5

  # x3 (active, uncorrelated) should have high inclusion
  expect_gt(incl[3], 0.7,
            label = paste("x3 inclusion:", round(incl[3], 3)))

  # x4, x5 (inactive) should have low inclusion
  expect_lt(incl[4], 0.3,
            label = paste("x4 inclusion:", round(incl[4], 3)))
  expect_lt(incl[5], 0.3,
            label = paste("x5 inclusion:", round(incl[5], 3)))

  # Key property: x1 + x2 combined mass captures the correlated signal.
  # Marginal x1 inclusion may be diluted by correlation with x2, so we
  # test the *pair* rather than x1 alone.
  expect_gt(incl[1] + incl[2], 0.8,
            label = paste("x1+x2 combined:", round(incl[1] + incl[2], 3)))

  # Model probabilities from discretized posterior draws surface the
  # structure that marginals obscure.
  draws <- posterior::as_draws_matrix(fit)
  b_cols <- sort(grep("^b_(?!Intercept)", colnames(draws), perl = TRUE, value = TRUE))
  b_draws <- as.matrix(draws[, b_cols])
  model_indicators <- apply(1L * (b_draws != 0), 1, paste0, collapse = "")
  model_freq <- sort(table(model_indicators), decreasing = TRUE)
  model_probs <- model_freq / sum(model_freq)

  # Column order after sort: b_x1, b_x2, b_x3, b_x4, b_x5
  # True model indicator: x1=on, x2=off, x3=on, x4=off, x5=off
  true_model <- "10100"
  # Competitor with swapped correlated pair: x1=off, x2=on, x3=on
  competitor <- "01100"

  # The true model should appear among the top models
  top_models <- names(model_probs)[1:min(5, length(model_probs))]
  expect_true(
    true_model %in% top_models || competitor %in% top_models,
    info = paste("Top 5 models:", paste(top_models, collapse = ", "),
                 "do not contain true or competitor model")
  )

  # Among models containing only one of {x1, x2}, the true model should
  # have higher or comparable probability to the competitor
  p_true <- if (true_model %in% names(model_probs)) model_probs[true_model] else 0
  p_comp <- if (competitor %in% names(model_probs)) model_probs[competitor] else 0
  expect_gte(as.numeric(p_true), as.numeric(p_comp) * 0.5,
             label = paste("P(true):", round(p_true, 4),
                           "P(competitor):", round(p_comp, 4)))
})

# ==============================================================================
# Section 9: Prior sensitivity (scale perturbation)
# ==============================================================================

test_that("brm_pdmp sticky: prior sensitivity — active > inactive under scale perturbation", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 100
  d <- 6
  X <- matrix(rnorm(n * d), n, d)
  colnames(X) <- paste0("x", 1:d)
  beta_true <- c(2.0, -1.5, rep(0, d - 2))
  y <- X %*% beta_true + rnorm(n, sd = 0.5)
  df <- data.frame(y = as.numeric(y), X)
  formula <- y ~ x1 + x2 + x3 + x4 + x5 + x6

  for (scale in c(2.5, 10)) {
    fit <- brm_pdmp(
      formula, data = df, family = gaussian(),
      prior = brms::set_prior(paste0("normal(0, ", scale, ")"), class = "b"),
      flow = "ZigZag", algorithm = "GridThinningStrategy",
      T = 100000, show_progress = FALSE,
      sticky = TRUE, model_prior = bernoulli(0.5)
    )
    incl <- attr(fit, "sticky")$inclusion_probs$chain1
    expect_true(min(incl[1:2]) > max(incl[3:d]),
                info = paste("Scale =", scale, "; active:", paste(round(incl[1:2], 3), collapse = ", "),
                             "; max inactive:", round(max(incl[3:d]), 3)))
  }
})
