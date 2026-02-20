# Tests for continuous-time estimators and parameter transforms.

skip_if_no_julia <- function() {
  julia_available <- tryCatch({
    JuliaCall::julia_setup(verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(julia_available, "Julia is not available")
}

# ─── S3 generic fallback tests (no Julia needed) ─────────────────────────────

test_that("var.default falls back to stats::var", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(var(x), stats::var(x))
})

test_that("sd.default falls back to stats::sd", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(sd(x), stats::sd(x))
})

test_that("cov.default falls back to stats::cov", {
  m <- matrix(rnorm(20), ncol = 2)
  expect_equal(cov(m), stats::cov(m))
})

test_that("cor.default falls back to stats::cor", {
  m <- matrix(rnorm(20), ncol = 2)
  expect_equal(cor(m), stats::cor(m))
})

# ─── Transform constructors (no Julia needed) ────────────────────────────────

test_that("identity_transform creates correct spec", {
  tr <- identity_transform()
  expect_equal(tr$type, "identity")
})

test_that("lower_transform creates correct spec", {
  tr <- lower_transform(0)
  expect_equal(tr$type, "lower")
  expect_equal(tr$lower, 0)
})

test_that("upper_transform creates correct spec", {
  tr <- upper_transform(10)
  expect_equal(tr$type, "upper")
  expect_equal(tr$upper, 10)
})

test_that("double_transform creates correct spec", {
  tr <- double_transform(0, 1)
  expect_equal(tr$type, "double")
  expect_equal(tr$lower, 0)
  expect_equal(tr$upper, 1)
})

test_that("double_transform rejects lower >= upper", {
  expect_error(double_transform(1, 0), "less than")
  expect_error(double_transform(1, 1), "less than")
})

test_that("transform constructors reject non-numeric input", {
  expect_error(lower_transform("a"), "type")
  expect_error(upper_transform("a"), "type")
  expect_error(double_transform("a", 1), "type")
})

# ─── pdmp_result class (no Julia needed) ─────────────────────────────────────

test_that("pdmp_result prints nicely", {
  r <- PDMPSamplersR:::new_pdmp_result(NULL, list(), 5L, 2L)
  expect_s3_class(r, "pdmp_result")
  expect_output(print(r), "5 dimensions.*2 chains")
})

test_that("estimator methods reject non-pdmp objects", {
  expect_error(mean.pdmp_result(42), "pdmp_result")
  expect_error(var(structure(list(), class = "wrong")), class = "simpleError")
})

# ─── Continuous-time estimator integration tests ─────────────────────────────

test_that("mean, var, sd produce correct-length output", {
  skip_if_no_julia()

  d <- 3
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 1000,
                        show_progress = FALSE)

  m <- mean(result)
  expect_length(m, d)
  expect_true(is.numeric(m))

  v <- var(result)
  expect_length(v, d)
  expect_true(all(v >= 0))

  s <- sd(result)
  expect_length(s, d)
  expect_true(all(s >= 0))
})

test_that("cov and cor produce correct-dimension output", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 1000,
                        show_progress = FALSE)

  C <- cov(result)
  expect_equal(dim(C), c(d, d))

  R <- cor(result)
  expect_equal(dim(R), c(d, d))
  # diagonal of correlation should be ~1
  expect_true(all(abs(diag(R) - 1) < 0.01))
})

test_that("quantile and median return correct types", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 1000,
                        show_progress = FALSE)

  q <- quantile(result, 0.5)
  expect_length(q, d)

  med <- median(result)
  expect_length(med, d)
})

test_that("ess returns positive values", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 2000,
                        show_progress = FALSE)

  e <- ess(result)
  expect_length(e, d)
  expect_true(all(e > 0))
})

test_that("cdf returns value in [0, 1]", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 1000,
                        show_progress = FALSE)

  p <- cdf(result, 0.0, coordinate = 1L)
  expect_true(is.numeric(p))
  expect_true(p >= 0 && p <= 1)
})

test_that("inclusion_probs returns correct-length output", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 1000,
                        show_progress = FALSE)

  ip <- inclusion_probs(result)
  expect_length(ip, d)
  # non-sticky sampler: coordinates are always nonzero, so inclusion ~ 1
  expect_true(all(is.finite(ip)))
})

# ─── Transformed estimators ─────────────────────────────────────────────────

test_that("mean with lower_transform produces positive values", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 2000,
                        show_progress = FALSE)

  transforms <- list(lower_transform(0), lower_transform(0))
  m <- mean(result, transforms = transforms)
  expect_length(m, d)
  # lower_transform(0) = exp(y), so mean should be > 0
  expect_true(all(m > 0))
})

test_that("mean with identity_transform matches unconstrained mean", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 2000,
                        show_progress = FALSE)

  transforms <- list(identity_transform(), identity_transform())
  m_identity <- mean(result, transforms = transforms)
  m_plain    <- mean(result)
  expect_equal(m_identity, m_plain, tolerance = 1e-10)
})

test_that("var with transforms returns non-negative values", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 2000,
                        show_progress = FALSE)

  transforms <- list(lower_transform(0), double_transform(0, 1))
  v <- var(result, transforms = transforms)
  expect_length(v, d)
  expect_true(all(v >= 0))
})

test_that("quantile with transforms respects monotonicity", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 2000,
                        show_progress = FALSE)

  transforms <- list(lower_transform(0), lower_transform(0))
  q25 <- quantile(result, 0.25, transforms = transforms, coordinate = 1L)
  q75 <- quantile(result, 0.75, transforms = transforms, coordinate = 1L)
  expect_true(q25 < q75)
  # lower_transform(0) output must be > 0
  expect_true(q25 > 0)
})

# ─── Discretize ─────────────────────────────────────────────────────────────

test_that("discretize with explicit dt works", {
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 500,
                        show_progress = FALSE)

  samples <- discretize(result, dt = 1.0)
  expect_true(is.matrix(samples))
  expect_equal(ncol(samples), d)
  # with dt = 1.0 and T = 500, expect ~500 rows
  expect_gt(nrow(samples), 100)
})
