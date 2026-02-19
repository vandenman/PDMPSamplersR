# Integration tests that require Julia to be set up.
# These tests are skipped if Julia is not available.

skip_if_no_julia <- function() {
  julia_available <- tryCatch({
    JuliaCall::julia_setup(verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(julia_available, "Julia is not available")

  # pdmp_available <- tryCatch({
  #   PDMPSamplersR:::check_for_julia_setup()
  #   # smoke test: run a trivial 1d ZigZag for 10 time units
  #   JuliaCall::julia_eval("let
  #     grad!(out, x) = (out .= x; out)
  #     r_pdmp_custom(grad!, 1, [0.0], \"ZigZag\", \"ThinningStrategy\", [0.0], ones(1,1);
  #                   t0=0.0, T=10.0, show_progress=false)
  #     true
  #   end")
  # }, error = function(e) FALSE)
  # testthat::skip_if_not(pdmp_available, "PDMPSamplers.jl Julia integration not functional")
}

test_that("pdmp_sample runs with a simple normal gradient", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 500,
                        show_progress = FALSE)

  expect_type(result, "list")
  samples <- result$samples
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

  expect_type(result, "list")
  samples <- result$samples
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

  expect_type(result, "list")
  samples <- result$samples
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

  expect_type(result, "list")
  expect_length(result, 2)
  for (chain in result) {
    expect_true(is.matrix(chain$samples))
    expect_equal(ncol(chain$samples), d)
  }
})

test_that("pdmp_sample posterior mean is reasonable for standard normal", {
  skip_on_cran()
  skip_if_no_julia()

  d <- 2
  neg_grad <- function(x) x

  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 5000,
                        show_progress = FALSE)

  samples <- result$samples
  posterior_means <- colMeans(samples)
  # for standard normal, posterior means should be near zero
  expect_true(all(abs(posterior_means) < 0.5))
})

test_that("pdmp_sample rejects invalid gradient function", {
  expect_error(pdmp_sample("not a function", d = 2, flow = "ZigZag", T = 100), "function")
  expect_error(pdmp_sample(function(x) x[1], d = 2, flow = "ZigZag", T = 100), "length")
})
