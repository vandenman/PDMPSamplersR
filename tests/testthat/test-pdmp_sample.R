# Integration tests that require Julia to be set up.
# These tests are skipped if Julia is not available.

test_that("pdmp_sample runs with a simple normal gradient", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 500,
                        show_progress = FALSE)

  expect_s3_class(result, "pdmp_result")
  expect_equal(result$d, d)
  expect_equal(result$n_chains, 1L)

  samples <- discretize(result)
  expect_true(is.matrix(samples))
  expect_equal(ncol(samples), d)
  expect_gt(nrow(samples), 10)
})

test_that("pdmp_sample with BouncyParticle flow", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 3
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "BouncyParticle", T = 500,
                        show_progress = FALSE)

  expect_s3_class(result, "pdmp_result")
  samples <- discretize(result)
  expect_true(is.matrix(samples))
  expect_equal(ncol(samples), d)
})

test_that("pdmp_sample with GridThinningStrategy", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x
  neg_hess <- function(x) diag(d)

  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag",
                        algorithm = "GridThinningStrategy", T = 500,
                        hessian = neg_hess, show_progress = FALSE)

  expect_s3_class(result, "pdmp_result")
  samples <- discretize(result)
  expect_true(is.matrix(samples))
  expect_equal(ncol(samples), d)
})

test_that("pdmp_sample with multiple chains", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 500,
                        n_chains = 2L, show_progress = FALSE)

  expect_s3_class(result, "pdmp_result")
  expect_equal(result$n_chains, 2L)

  samples_1 <- discretize(result, chain = 1L)
  samples_2 <- discretize(result, chain = 2L)
  expect_true(is.matrix(samples_1))
  expect_true(is.matrix(samples_2))
  expect_equal(ncol(samples_1), d)
  expect_equal(ncol(samples_2), d)
})

test_that("pdmp_sample posterior mean is reasonable for standard normal", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 5000,
                        show_progress = FALSE)

  posterior_means <- mean(result)
  expect_length(posterior_means, d)
  expect_true(all(abs(posterior_means) < 0.5))
})

test_that("pdmp_sample rejects invalid gradient function", {
  expect_error(pdmp_sample("not a function", d = 2, flow = "ZigZag", T = 100), "function")
  expect_error(pdmp_sample(function(x) x[1], d = 2, flow = "ZigZag", T = 100), "length")
})

test_that("pdmp_sample forwards support-boundary diagnostics", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) {
    if (x[1] >= 1) stop("Outside support")
    x
  }
  neg_hess <- function(x) diag(d)

  expect_error(
    pdmp_sample(
      neg_grad, d = d, flow = "BouncyParticle",
      algorithm = "GridThinningStrategy", T = 10,
      x0 = c(0, 0), theta0 = c(1, 0), hessian = neg_hess,
      show_progress = FALSE, materialize = FALSE,
      support_boundary = support_boundary_control(
        mode = "line_search",
        max_bisection_steps = 3
      )
    ),
    "SupportBoundaryError|Boundary localized",
    class = "pdmp_support_boundary_error"
  )
})

test_that("pdmp_sample accepts line_search_truncated_refresh for BPS-family flows", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) {
    if (x[1] >= 1) stop("Outside support")
    x
  }
  neg_hess <- function(x) diag(d)

  result <- pdmp_sample(
    neg_grad, d = d, flow = "BouncyParticle",
    algorithm = "GridThinningStrategy", T = 5,
    x0 = c(0, 0), theta0 = c(1, 0), hessian = neg_hess,
    show_progress = FALSE, materialize = FALSE,
    support_boundary = support_boundary_control(
      mode = "line_search_truncated_refresh",
      max_refresh_attempts = 200,
      refresh_probe_time = 1e-4
    )
  )

  expect_s3_class(result, "pdmp_result")
  expect_gte(result$stats$support_boundary_events[[1]], 1)
})

# ──────────────────────────────────────────────────────────────────────────────
# Smoke tests for all dynamics (just check they run, not correctness)
# ──────────────────────────────────────────────────────────────────────────────

test_that("pdmp_sample runs with Boomerang flow", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "Boomerang",
                        algorithm = "GridThinningStrategy", T = 500,
                        show_progress = FALSE)

  expect_s3_class(result, "pdmp_result")
  expect_equal(result$d, d)
})

test_that("pdmp_sample runs with AdaptiveBoomerang (diagonal)", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "AdaptiveBoomerang",
                        algorithm = "GridThinningStrategy", T = 5000,
                        adaptive_scheme = "diagonal",
                        show_progress = FALSE)

  expect_s3_class(result, "pdmp_result")
  expect_equal(result$d, d)
})

test_that("pdmp_sample runs with AdaptiveBoomerang (fullrank)", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "AdaptiveBoomerang",
                        algorithm = "GridThinningStrategy", T = 5000,
                        adaptive_scheme = "fullrank",
                        show_progress = FALSE)

  expect_s3_class(result, "pdmp_result")
  expect_equal(result$d, d)
})

test_that("pdmp_sample runs with PreconditionedZigZag", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "PreconditionedZigZag",
                        algorithm = "GridThinningStrategy", T = 500,
                        show_progress = FALSE)

  expect_s3_class(result, "pdmp_result")
  expect_equal(result$d, d)
})

test_that("pdmp_sample runs with PreconditionedBPS", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "PreconditionedBPS",
                        algorithm = "GridThinningStrategy", T = 500,
                        show_progress = FALSE)

  expect_s3_class(result, "pdmp_result")
  expect_equal(result$d, d)
})
