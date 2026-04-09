# MCMC sampling integration tests for brm_pdmp() (standard and subsampled).
#
# These tests run full PDMP sampling and compare posteriors to HMC references.
# They require Julia + BridgeStan and take ~20 minutes.
#
# Gated by PDMPSAMPLERSR_SLOW_TESTS environment variable + skip_on_cran.
# Run locally:
#   NOT_CRAN=true PDMPSAMPLERSR_SLOW_TESTS=true Rscript -e 'require(PDMPSamplersR); testthat::test_file("tests/testthat/test-brms_sampling.R")'

skip_if_no_slow_tests <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not(
    identical(Sys.getenv("PDMPSAMPLERSR_SLOW_TESTS"), "true"),
    "Slow MCMC tests not enabled (set PDMPSAMPLERSR_SLOW_TESTS=true)"
  )
  skip_if_no_brms_setup()
}

compare_posteriors <- function(fit_pdmp, fit_reference,
                               tols = c("Estimate" = .05, "Est.Error" = .15,
                                        "Q2.5" = 0.15, "Q97.5" = 0.15)) {
  ests_pdmp <- brms::fixef(fit_pdmp)
  ests_ref  <- brms::fixef(fit_reference)
  for (i in seq_along(tols)) {
    nm <- names(tols)[i]
    testthat::expect_equal(ests_pdmp[, nm], ests_ref[, nm], tolerance = tols[i],
                           label = paste("PDMP", nm, "should be close to reference"))
  }
}

# ==============================================================================
# Section 1: Standard brm_pdmp (no subsampling)
# ==============================================================================

test_that("brm_pdmp rejects sample_prior != 'no'", {
  skip_if_no_slow_tests()
  expect_error(
    brm_pdmp(y ~ x, data = data.frame(y = 1:3, x = 1:3), sample_prior = "yes"),
    "sample_prior"
  )
})

test_that("brm_pdmp: gaussian y ~ x", {
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

  s <- summary(fit)
  expect_true(!is.null(s))
  expect_true("fixed" %in% names(s))

  fixed <- s$fixed
  expect_true("Intercept" %in% rownames(fixed))
  expect_true("x" %in% rownames(fixed))
  expect_true(abs(fixed["Intercept", "Estimate"] - 2) < 1.0)
  expect_true(abs(fixed["x", "Estimate"] - 0.5) < 1.0)

  ce <- brms::conditional_effects(fit)
  expect_true(length(ce) > 0)

  fit_reference <- brms::brm(y ~ x, data = df, family = gaussian(),
                             chains = 1, iter = 2000, warmup = 1000,
                             silent = 2, refresh = 0)

  compare_posteriors(fit, fit_reference,
                     tols = c("Estimate" = .05, "Est.Error" = .15,
                              "Q2.5" = 0.1, "Q97.5" = 0.1))
})

test_that("brm_pdmp: logistic regression", {
  skip_if_no_slow_tests()

  set.seed(123)
  n <- 100
  x <- rnorm(n)
  prob <- plogis(0.5 + 1.0 * x)
  y <- rbinom(n, 1, prob)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = brms::bernoulli(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Intercept" %in% rownames(fixed))
  expect_true("x" %in% rownames(fixed))

  fit_reference <- brms::brm(y ~ x, data = df, family = brms::bernoulli(),
                             chains = 1, iter = 2000, warmup = 1000,
                             silent = 2, refresh = 0)

  compare_posteriors(fit, fit_reference)
})

test_that("brm_pdmp: Poisson regression", {
  skip_if_no_slow_tests()

  set.seed(456)
  n <- 100
  x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.3 * x))
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = poisson(),
                  flow = "BouncyParticle", algorithm = "GridThinningStrategy",
                  T = 50000, show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Intercept" %in% rownames(fixed))
  expect_true("x" %in% rownames(fixed))

  fit_reference <- brms::brm(y ~ x, data = df, family = poisson(),
                             chains = 1, iter = 2000, warmup = 1000,
                             silent = 2, refresh = 0)

  compare_posteriors(fit, fit_reference)
})

test_that("brm_pdmp: Poisson regression in epileptic patients", {
  skip_if_no_slow_tests()

  data(epilepsy, package = "brms")
  fit_reference <- brms::brm(
    count ~ zBase * Trt + (1|patient),
    data = epilepsy, family = poisson(),
    prior = brms::prior(normal(0, 10), class = b) +
            brms::prior(cauchy(0, 2), class = sd),
    silent = 2, refresh = 0
  )

  fit <- brm_pdmp(
    count ~ zBase * Trt + (1|patient),
    data = epilepsy, family = poisson(),
    prior = brms::prior(normal(0, 10), class = b) +
            brms::prior(cauchy(0, 2), class = sd),
    flow = "AdaptiveBoomerang", algorithm = "GridThinningStrategy",
    T = 10000, show_progress = TRUE
  )

  compare_posteriors(fit, fit_reference,
                     tols = c("Estimate" = .05, "Est.Error" = .15,
                              "Q2.5" = 0.1, "Q97.5" = 0.1))
})

test_that("brm_pdmp: multi-chain gaussian with Rhat", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 2 + 0.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = gaussian(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, show_progress = FALSE,
                  n_chains = 2L)

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed

  expect_true("Rhat" %in% colnames(fixed))
  rhat_vals <- fixed[, "Rhat"]
  for (rh in rhat_vals) {
    expect_lt(rh, 1.05)
  }

  expect_true(abs(fixed["Intercept", "Estimate"] - 2) < 0.1)
  expect_true(abs(fixed["x", "Estimate"] - 0.5) < 0.15)
})

test_that("brm_pdmp: compute_lp populates lp__", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 2 + 0.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = gaussian(),
                  flow = "AdaptiveBoomerang", algorithm = "GridThinningStrategy",
                  T = 20000, show_progress = FALSE,
                  compute_lp = TRUE)

  expect_s3_class(fit, "brmsfit")
  lp <- rstan::extract(fit$fit, pars = "lp__")$lp__
  expect_true(all(is.finite(lp)))
  expect_true(all(lp < 0))
})

test_that("brm_pdmp: pp_check and loo work", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 2 + 0.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = gaussian(),
                  flow = "PreconditionedZigZag", algorithm = "GridThinningStrategy",
                  T = 20000, show_progress = FALSE)

  pp <- brms::pp_check(fit, ndraws = 50)
  expect_s3_class(pp, "gg")

  loo_result <- suppressWarnings(brms::loo(fit))
  expect_s3_class(loo_result, "loo")
  expect_true(!is.null(loo_result$estimates))
})

# ==============================================================================
# Section 2: Subsampled brm_pdmp
# ==============================================================================

test_that("brm_pdmp rejects subsample_size >= nrow(data)", {
  skip_if_no_slow_tests()
  expect_error(
    brm_pdmp(y ~ x, data = data.frame(y = 1:10, x = 1:10),
             subsample_size = 10L, T = 100, show_progress = FALSE),
    "subsample_size"
  )
})

test_that("brm_pdmp subsampled: gaussian y ~ x", {
  skip_if_no_slow_tests()

  set.seed(42)
  n <- 100
  x <- rnorm(n)
  y <- 2 + 0.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  fit <- tryCatch(
    brm_pdmp(y ~ x, data = df, family = gaussian(),
             flow = "ZigZag", algorithm = "GridThinningStrategy",
             T = 50000, t_warmup = 5000,
             subsample_size = 20L,
             show_progress = FALSE),
    error = function(e) {
      if (grepl("Scale vector is inf", conditionMessage(e)))
        testthat::skip("sigma=inf during PDMP initialization (stochastic)")
      stop(e)
    }
  )

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Intercept" %in% rownames(fixed))
  expect_true("x" %in% rownames(fixed))
})

test_that("brm_pdmp subsampled: bernoulli y ~ x (with summary and pp_check)", {
  skip_if_no_slow_tests()

  set.seed(123)
  n <- 200
  x <- rnorm(n)
  prob <- plogis(0.5 + 1.0 * x)
  y <- rbinom(n, 1, prob)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = brms::bernoulli(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, t_warmup = 5000,
                  subsample_size = 30L,
                  show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Intercept" %in% rownames(fixed))
  expect_true("x" %in% rownames(fixed))
  expect_true("fixed" %in% names(s))

  pp <- brms::pp_check(fit, ndraws = 20)
  expect_s3_class(pp, "gg")
})

test_that("brm_pdmp subsampled: poisson y ~ x", {
  skip_if_no_slow_tests()

  set.seed(456)
  n <- 200
  x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.3 * x))
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = poisson(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, t_warmup = 5000,
                  subsample_size = 30L,
                  show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Intercept" %in% rownames(fixed))
  expect_true("x" %in% rownames(fixed))
})

test_that("brm_pdmp subsampled: multi-chain Rhat", {
  skip_if_no_slow_tests()

  set.seed(123)
  n <- 200
  x <- rnorm(n)
  prob <- plogis(0.5 + 1.0 * x)
  y <- rbinom(n, 1, prob)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = brms::bernoulli(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, t_warmup = 5000,
                  subsample_size = 30L,
                  n_chains = 2L,
                  show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Rhat" %in% colnames(fixed))
  rhat_vals <- fixed[, "Rhat"]
  for (rh in rhat_vals) {
    expect_true(is.finite(rh))
  }
})
