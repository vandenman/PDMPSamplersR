# Integration tests that require Julia to be set up.
# These tests are skipped if Julia is not available.

if (!exists("with_seed", envir = asNamespace("PDMPSamplersR"), inherits = FALSE) &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
}

pdmp_sample <- PDMPSamplersR::pdmp_sample
with_seed <- get("with_seed", envir = asNamespace("PDMPSamplersR"), inherits = FALSE)

expect_identical_skeleton <- function(result_1, result_2) {
  expect_identical(result_1$skeleton[[1]]$times, result_2$skeleton[[1]]$times)
  expect_identical(result_1$skeleton[[1]]$positions, result_2$skeleton[[1]]$positions)
  expect_identical(result_1$skeleton[[1]]$velocities, result_2$skeleton[[1]]$velocities)
}

test_that("with_seed evaluates code reproducibly and restores R RNG state", {
  set.seed(2026)
  seed_before <- .Random.seed
  kind_before <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(kind_before))
    assign(".Random.seed", seed_before, envir = .GlobalEnv)
  }, add = TRUE)

  draw_1 <- with_seed(stats::rnorm(5), seed = 7)
  seed_after <- .Random.seed
  kind_after <- RNGkind()
  draw_2 <- with_seed(stats::rnorm(5), seed = 7)

  expect_identical(draw_1, draw_2)
  expect_identical(seed_after, seed_before)
  expect_identical(kind_after, kind_before)

  rm(".Random.seed", envir = .GlobalEnv)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  invisible(with_seed(stats::rnorm(1), seed = 8))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("pdmp_sample with a seed returns an identical ZigZag trace twice", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  d <- 2
  neg_grad <- function(x) x
  x0 <- c(0.25, -0.75)

  result_1 <- pdmp_sample(
    neg_grad, d = d, flow = "ZigZag", T = 500,
    x0 = x0, seed = 42, show_progress = FALSE
  )
  result_2 <- pdmp_sample(
    neg_grad, d = d, flow = "ZigZag", T = 500,
    x0 = x0, seed = 42, show_progress = FALSE
  )

  expect_identical_skeleton(result_1, result_2)
})

test_that("pdmp_sample with a seed returns an identical AdaptiveBoomerang trace twice", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  d <- 2
  neg_grad <- function(x) x

  result_1 <- pdmp_sample(
    neg_grad, d = d,
    flow = "AdaptiveBoomerang",
    T = 5000,
    algorithm = "GridThinningStrategy",
    adaptive_scheme = "fullrank",
    seed = 1234,
    show_progress = FALSE
  )

  result_2 <- pdmp_sample(
    neg_grad, d = d,
    flow = "AdaptiveBoomerang",
    T = 5000,
    algorithm = "GridThinningStrategy",
    adaptive_scheme = "fullrank",
    seed = 1234,
    show_progress = FALSE
  )

  expect_identical_skeleton(result_1, result_2)
})

test_that("pdmp_sample seed controls default x0 without advancing R's RNG", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  d <- 2
  neg_grad <- function(x) x

  set.seed(2026)
  seed_before <- .Random.seed

  result_1 <- pdmp_sample(
    neg_grad, d = d, flow = "ZigZag", T = 500,
    seed = 7, show_progress = FALSE
  )
  seed_after <- .Random.seed

  result_2 <- pdmp_sample(
    neg_grad, d = d, flow = "ZigZag", T = 500,
    seed = 7, show_progress = FALSE
  )

  expect_identical(seed_after, seed_before)
  expect_identical_skeleton(result_1, result_2)
})
