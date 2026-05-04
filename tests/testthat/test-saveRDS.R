# Tests for saveRDS / materialize round-trip support.

# ─── materialize error without Julia (no Julia needed) ───────────────────────

test_that("estimator on bare save (no skeleton) gives informative error", {
  result <- PDMPSamplersR:::new_pdmp_result(NULL, list(), 2L, 1L)
  expect_error(mean(result), "skeleton")
})

# ─── Round-trip tests (Julia required) ───────────────────────────────────────

.run_roundtrip <- function(flow, d = 2, T = 2000) {
  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = d, flow = flow, T = T,
                        show_progress = FALSE)
  list(
    original = result,
    reloaded = {
      m <- materialize(result)
      tmp <- tempfile(fileext = ".rds")
      on.exit(unlink(tmp))
      saveRDS(m, tmp)
      readRDS(tmp)
    }
  )
}

test_that("ZigZag round-trip preserves estimators", {
  skip_if_no_pdmp_julia_backend()

  pair <- .run_roundtrip("ZigZag")
  orig <- pair$original
  rel  <- pair$reloaded

  expect_equal(mean(rel),                  mean(orig),                  tolerance = 1e-8)
  expect_equal(var(rel),                   var(orig),                   tolerance = 1e-8)
  expect_equal(cov(rel),                   cov(orig),                   tolerance = 1e-8)
  expect_equal(ess(rel),                   ess(orig),                   tolerance = 1e-8)
  expect_equal(cdf(rel, 0, coordinate = 1L), cdf(orig, 0, coordinate = 1L), tolerance = 1e-8)
  expect_equal(quantile(rel, 0.5),         quantile(orig, 0.5),         tolerance = 1e-8)
  expect_equal(discretize(rel),            discretize(orig),            tolerance = 1e-8)
})

test_that("Boomerang round-trip preserves estimators", {
  skip_if_no_pdmp_julia_backend()

  pair <- .run_roundtrip("Boomerang")
  orig <- pair$original
  rel  <- pair$reloaded

  expect_equal(mean(rel),                   mean(orig),                   tolerance = 1e-8)
  expect_equal(cdf(rel, 0, coordinate = 1L), cdf(orig, 0, coordinate = 1L), tolerance = 1e-8)
})

test_that("AdaptiveBoomerang round-trip preserves estimators", {
  skip_if_no_pdmp_julia_backend()

  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = 2, flow = "AdaptiveBoomerang",
                        algorithm = "GridThinningStrategy",
                        T = 3000, t_warmup = 500, show_progress = FALSE)

  m <- materialize(result)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(m, tmp)
  rel <- readRDS(tmp)

  expect_equal(mean(rel),                   mean(result),                   tolerance = 1e-8)
  expect_equal(cdf(rel, 0, coordinate = 1L), cdf(result, 0, coordinate = 1L), tolerance = 1e-8)
})

test_that("transformed estimators match after ZigZag round-trip", {
  skip_if_no_pdmp_julia_backend()

  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = 2, flow = "ZigZag", T = 2000,
                        show_progress = FALSE)
  tr <- list(lower_transform(0), lower_transform(0))

  m <- materialize(result)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(m, tmp)
  rel <- readRDS(tmp)

  expect_equal(mean(rel, transforms = tr), mean(result, transforms = tr), tolerance = 1e-8)
})

test_that("materialize is idempotent", {
  skip_if_no_pdmp_julia_backend()

  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = 2, flow = "ZigZag", T = 1000,
                        show_progress = FALSE)

  m1 <- materialize(result)
  m2 <- materialize(m1)

  expect_equal(m1$skeleton, m2$skeleton)
  expect_null(m2$chains)
})
