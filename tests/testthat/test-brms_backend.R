# Integration tests for brms backend.
# Skipped if brms, rstan, or Julia are not available.

skip_if_no_brms_setup <- function() {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("rstan")

  julia_available <- tryCatch({
    PDMPSamplersR:::check_for_julia_setup()
    JuliaCall::julia_eval("hasmethod(PDMPModel, Tuple{String, String})")
  }, error = function(e) FALSE)
  testthat::skip_if_not(julia_available, "Julia + BridgeStan not available")
}

test_that("brm_pdmp rejects sample_prior != 'no'", {
  skip_if_no_brms_setup()
  expect_error(
    brm_pdmp(y ~ x, data = data.frame(y = 1:3, x = 1:3), sample_prior = "yes"),
    "sample_prior"
  )
})

test_that("brm_pdmp: gaussian y ~ x", {
  skip_if_no_brms_setup()

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

  # Intercept should be near 2
  expect_true(abs(fixed["Intercept", "Estimate"] - 2) < 1.0)
  # Slope should be near 0.5
  expect_true(abs(fixed["x", "Estimate"] - 0.5) < 1.0)

  # conditional_effects should work
  ce <- brms::conditional_effects(fit)
  expect_true(length(ce) > 0)

  fit_reference <- brms::brm(y ~ x, data = df, family = gaussian(),
                             chains = 1, iter = 2000, warmup = 1000, silent = 2, refresh = 0)

  ests_reference <- brms::fixef(fit_reference)
  ests_pdmp      <- brms::fixef(fit)

  tols <- c("Estimate" = .05, "Est.Error"  = .15, "Q2.5" = 0.1, "Q97.5" = 0.1)
  for (i in seq_along(tols)) {
    nm <- names(tols)[i]
    testthat::expect_equal(ests_pdmp[, nm], ests_reference[, nm], tolerance = tols[i],
                           label = paste("PDMP", nm, "should be close to reference"))
  }
})

test_that("brm_pdmp: logistic regression", {
  skip_if_no_brms_setup()

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
                             chains = 1, iter = 2000, warmup = 1000, silent = 2, refresh = 0)

  ests_reference <- brms::fixef(fit_reference)
  ests_pdmp      <- brms::fixef(fit)

  tols <- c("Estimate" = .05, "Est.Error"  = .15, "Q2.5" = 0.1, "Q97.5" = 0.1)
  for (i in seq_along(tols)) {
    nm <- names(tols)[i]
    testthat::expect_equal(ests_pdmp[, nm], ests_reference[, nm], tolerance = tols[i],
                           label = paste("PDMP", nm, "should be close to reference"))
  }
})

test_that("brm_pdmp: Poisson regression", {
  skip_if_no_brms_setup()

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
                             chains = 1, iter = 2000, warmup = 1000, silent = 2, refresh = 0)

  ests_reference <- brms::fixef(fit_reference)
  ests_pdmp      <- brms::fixef(fit)

  tols <- c("Estimate" = .05, "Est.Error"  = .15, "Q2.5" = 0.1, "Q97.5" = 0.1)
  for (i in seq_along(tols)) {
    nm <- names(tols)[i]
    testthat::expect_equal(ests_pdmp[, nm], ests_reference[, nm], tolerance = tols[i],
                           label = paste("PDMP", nm, "should be close to reference"))
  }
})

test_that("brm_pdmp: Poisson regression example in epileptic patients", {
  skip_if_no_brms_setup()

  data(epilepsy, package = "brms")
  fit_reference <- brms::brm(
    count ~ zBase * Trt + (1|patient),
    data = epilepsy, family = poisson(),
    prior = brms::prior(normal(0, 10), class = b) + brms::prior(cauchy(0, 2), class = sd),
    silent = 2, refresh = 0
  )

  fit <- brm_pdmp(
    count ~ zBase * Trt + (1|patient),
    data = epilepsy, family = poisson(),
    prior = brms::prior(normal(0, 10), class = b) + brms::prior(cauchy(0, 2), class = sd),
    flow = "AdaptiveBoomerang", algorithm = "GridThinningStrategy",
    T = 10000, show_progress = TRUE
  )

  # generate a summary of the results
  # summary(fit_reference)
  # summary(fit)

  # plot the MCMC chains as well as the posterior distributions
  # plot(fit_reference)
  # plot(fit)

  ests_reference <- brms::fixef(fit_reference)
  ests_pdmp      <- brms::fixef(fit)

  tols <- c("Estimate" = .05, "Est.Error"  = .15, "Q2.5" = 0.1, "Q97.5" = 0.1)
  for (i in seq_along(tols)) {
    nm <- names(tols)[i]
    testthat::expect_equal(ests_pdmp[, nm], ests_reference[, nm], tolerance = tols[i],
                           label = paste("PDMP", nm, "should be close to reference"))
  }
})

