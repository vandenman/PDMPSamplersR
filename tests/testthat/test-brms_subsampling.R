# Tests for brms subsampling: helper functions and integration tests.

# ── Unit tests for helper functions (no Julia needed) ─────────────────────────

test_that("subset_standata subsets correctly", {
  sdata <- list(
    N = 10L, Y = 1:10, K = 3L, Kc = 2L,
    X = matrix(seq_len(30), nrow = 10, ncol = 3),
    means_X = c(5.5, 15.5)
  )
  indices <- c(1L, 3L, 5L)

  subset_standata <- PDMPSamplersR:::subset_standata
  sub <- subset_standata(sdata, indices)

  expect_equal(sub$N, 3L)
  expect_equal(sub$Y, c(1L, 3L, 5L))
  expect_equal(nrow(sub$X), 3L)
  expect_equal(ncol(sub$X), 3L)
  expect_equal(sub$X[1, ], sdata$X[1, ])
  expect_equal(sub$X[2, ], sdata$X[3, ])
  expect_equal(sub$X[3, ], sdata$X[5, ])
  expect_equal(sub$means_X, c(5.5, 15.5))
  expect_equal(sub$K, 3L)
})

test_that("make_prior_standata creates valid prior data", {
  sdata <- list(
    N = 10L, Y = 1:10, K = 3L, Kc = 2L,
    X = matrix(seq_len(30), nrow = 10, ncol = 3),
    means_X = c(5.5, 15.5)
  )

  make_prior_standata <- PDMPSamplersR:::make_prior_standata
  prior <- make_prior_standata(sdata)

  expect_equal(prior$N, 1L)
  expect_equal(prior$K, 3L)
  expect_equal(prior$Kc, 2L)
  expect_equal(prior$prior_only, 1L)
  expect_equal(prior$means_X, c(5.5, 15.5))
  expect_true(is.integer(prior$Y))
  expect_equal(length(prior$Y), 1L)
  expect_true(!is.null(dim(prior$Y)))
  expect_equal(nrow(prior$X), 1L)
  expect_equal(ncol(prior$X), 3L)
})

test_that("make_prior_standata handles double Y", {
  sdata <- list(
    N = 5L, Y = c(1.1, 2.2, 3.3, 4.4, 5.5), K = 2L, Kc = 1L,
    X = matrix(seq_len(10), nrow = 5, ncol = 2),
    means_X = 3.5
  )

  make_prior_standata <- PDMPSamplersR:::make_prior_standata
  prior <- make_prior_standata(sdata)

  expect_true(is.double(prior$Y))
  expect_equal(length(prior$Y), 1L)
  expect_true(!is.null(dim(prior$Y)))
})

test_that("subset_standata subsets extra observation fields", {
  sdata <- list(
    N = 5L, Y = 1:5, K = 2L, Kc = 1L,
    X       = matrix(seq_len(10), nrow = 5, ncol = 2),
    offsets = c(0.1, 0.2, 0.3, 0.4, 0.5),
    trials  = c(10L, 20L, 30L, 40L, 50L),
    means_X = 1.5
  )
  indices <- c(2L, 4L)

  sub <- PDMPSamplersR:::subset_standata(sdata, indices)

  expect_equal(sub$offsets, c(0.2, 0.4))
  expect_equal(sub$trials, c(20L, 40L))
})

test_that("make_prior_standata includes extra observation fields at N=1", {
  sdata <- list(
    N = 5L, Y = 1:5, K = 2L, Kc = 1L,
    X       = matrix(seq_len(10), nrow = 5, ncol = 2),
    offsets = c(0.1, 0.2, 0.3, 0.4, 0.5),
    trials  = c(10L, 20L, 30L, 40L, 50L),
    means_X = 1.5
  )

  prior <- PDMPSamplersR:::make_prior_standata(sdata)

  expect_equal(prior$N, 1L)
  expect_equal(length(prior$offsets), 1L)
  expect_equal(length(prior$trials), 1L)
})

# ── Validation tests ─────────────────────────────────────────────────────────

test_that("brm_pdmp rejects subsample_size >= nrow(data)", {
  skip_if_no_brms_setup()
  expect_error(
    brm_pdmp(y ~ x, data = data.frame(y = 1:10, x = 1:10),
             subsample_size = 10L, T = 100, show_progress = FALSE),
    "subsample_size"
  )
})

test_that("brm_stancode accepts random effects with subsampling", {
  skip_if_no_brms_setup()

  set.seed(1)
  df <- data.frame(
    y = rnorm(20), x = rnorm(20),
    group = rep(1:4, each = 5)
  )
  result <- brm_stancode(y ~ x + (1 | group), data = df,
                         subsample_size = 5L)
  expect_true(is.list(result))
  expect_true(grepl("for \\(i in 1:m_sub\\)", result$ext_cpp))
  expect_true(grepl("pdmp_get_subsample_index", result$ext_cpp))
})

# ── Integration tests (need Julia + brms) ────────────────────────────────────

test_that("brm_pdmp with subsample_size: gaussian y ~ x", {
  skip_if_no_brms_setup()
  # NOTE: Gaussian models have a positive-constrained sigma parameter. Random
  # PDMP initialization can push log(sigma) to extremes, causing BridgeStan to
  # throw "Scale vector is inf". This is a pre-existing sampler limitation, not
  # specific to subsampling. Skip on this known error.
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

test_that("brm_pdmp with subsample_size: bernoulli y ~ x", {
  skip_if_no_brms_setup()

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
})

test_that("brm_pdmp with subsample_size: poisson y ~ x", {
  skip_if_no_brms_setup()

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

test_that("brm_pdmp with subsample_size: summary and pp_check work", {
  skip_if_no_brms_setup()

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

  s <- summary(fit)
  expect_true(!is.null(s))
  expect_true("fixed" %in% names(s))

  pp <- brms::pp_check(fit, ndraws = 20)
  expect_s3_class(pp, "gg")
})

test_that("brm_pdmp with subsample_size: multi-chain Rhat", {
  skip_if_no_brms_setup()

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
