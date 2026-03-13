# Integration tests for brms backend.
# Skipped if brms, rstan, or Julia are not available.

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

  tols <- c("Estimate" = .05, "Est.Error"  = .15, "Q2.5" = 0.15, "Q97.5" = 0.1)
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

  tols <- c("Estimate" = .05, "Est.Error"  = .15, "Q2.5" = 0.15, "Q97.5" = 0.15)
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


test_that("brm_pdmp: multi-chain gaussian with Rhat", {
  skip_if_no_brms_setup()

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
  skip_if_no_brms_setup()

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

test_that("brm_pdmp: student_t family", {
  skip_if_no_brms_setup()

  set.seed(100)
  n <- 80
  x <- rnorm(n)
  y <- 1.5 + 0.8 * x + rt(n, df = 3)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = brms::student(),
                  flow = "AdaptiveBoomerang", algorithm = "GridThinningStrategy",
                  T = 50000, show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")
  s <- summary(fit)
  expect_true("Intercept" %in% rownames(s$fixed))
  expect_true(abs(s$fixed["Intercept", "Estimate"] - 1.5) < 0.1)
  expect_true(abs(s$fixed["x", "Estimate"] - 0.8) < 0.15)
})

test_that("brm_pdmp: negbinomial family", {
  skip_if_no_brms_setup()

  set.seed(200)
  n <- 100
  x <- rnorm(n)
  mu <- exp(0.5 + 0.3 * x)
  y <- rnbinom(n, size = 5, mu = mu)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = brms::negbinomial(),
                  flow = "AdaptiveBoomerang", algorithm = "GridThinningStrategy",
                  T = 50000, show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")
  s <- summary(fit)
  expect_true("Intercept" %in% rownames(s$fixed))
  expect_true(abs(s$fixed["Intercept", "Estimate"] - 0.5) < 1.0)
  expect_true(abs(s$fixed["x", "Estimate"] - 0.3) < 1.0)
})

test_that("brm_pdmp: Beta family", {
  skip_if_no_brms_setup()

  set.seed(300)
  n <- 100
  x <- rnorm(n)
  mu <- plogis(0.5 + 0.4 * x)
  phi <- 10
  y <- rbeta(n, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  y <- pmin(pmax(y, 0.001), 0.999) # avoid (unlikely) exact 0 or 1 values which cause issues
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = brms::Beta(),
                  flow = "AdaptiveBoomerang", algorithm = "GridThinningStrategy",
                  T = 50000, show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")
  s <- summary(fit)
  expect_true("Intercept" %in% rownames(s$fixed))
  expect_true("x" %in% rownames(s$fixed))
})

test_that("brm_pdmp: pp_check and loo work", {
  skip_if_no_brms_setup()

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
  # TODO: some substantial check on loo result?
})
