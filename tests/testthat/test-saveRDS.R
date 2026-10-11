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

  neg_grad <- function(x) c(x[[1L]] / 4, x[[2L]] * 2)
  result <- pdmp_sample(neg_grad, d = 2, flow = "AdaptiveBoomerang",
                        algorithm = "GridThinningStrategy",
                        T = 3000, t_warmup = 500, show_progress = FALSE)

  m <- materialize(result)
  expect_equal(dim(m$skeleton[[1L]]$flow_cov), c(2L, 2L))
  expect_true(isSymmetric(m$skeleton[[1L]]$flow_cov))
  expect_true(all(eigen(m$skeleton[[1L]]$flow_cov, symmetric = TRUE)$values > 0))
  expect_false(isTRUE(all.equal(m$skeleton[[1L]]$flow_cov, diag(2))))
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
  expect_identical(m1$initial_states, m2$initial_states)
  expect_identical(m1$retained_initial_states, m2$retained_initial_states)
  expect_identical(m1$terminal_states, m2$terminal_states)
  expect_null(m2$chains)
})

# ─── Sparse skeleton (ZigZag) ─────────────────────────────────────────────────

test_that("ZigZag materialize produces a sparse skeleton", {
  skip_if_no_pdmp_julia_backend()

  neg_grad <- function(x) x
  result <- pdmp_sample(neg_grad, d = 2, flow = "ZigZag", T = 1000,
                        show_progress = FALSE)
  m <- materialize(result)

  expect_true(m$skeleton[[1]]$sparse)
  expect_true(!is.null(m$skeleton[[1]]$event_indices))
  expect_true(is.null(m$skeleton[[1]]$times))
})

test_that("ZigZag sparse round-trip preserves estimators and is O(N) not O(d x N)", {
  skip_if_no_pdmp_julia_backend()

  neg_grad <- function(x) x
  d <- 4L
  result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 3000,
                        show_progress = FALSE)

  m <- materialize(result)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(m, tmp)
  rel <- readRDS(tmp)

  expect_equal(mean(rel),                    mean(result),                    tolerance = 1e-8)
  expect_equal(var(rel),                     var(result),                     tolerance = 1e-8)
  expect_equal(cov(rel),                     cov(result),                     tolerance = 1e-8)
  expect_equal(ess(rel),                     ess(result),                     tolerance = 1e-8)
  expect_equal(cdf(rel, 0, coordinate = 1L), cdf(result, 0, coordinate = 1L), tolerance = 1e-8)

  n_events <- length(m$skeleton[[1]]$event_times)
  sparse_size  <- 4L * n_events           # index + time + position + velocity per event
  dense_size   <- (2L * d + 1L) * n_events  # d positions + d velocities + 1 time per event
  expect_lt(sparse_size, dense_size)
})

test_that("Boomerang with non-zero centering round-trip preserves estimators", {
  skip_if_no_pdmp_julia_backend()

  mu <- c(2.0, -1.0)
  neg_grad <- function(x) x - mu
  result <- pdmp_sample(neg_grad, d = 2, flow = "Boomerang", T = 2000,
                        flow_mean = mu, show_progress = FALSE)

  m <- materialize(result)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(m, tmp)
  rel <- readRDS(tmp)

  expect_equal(mean(rel),                    mean(result),                    tolerance = 1e-8)
  expect_equal(var(rel),                     var(result),                     tolerance = 1e-8)
  expect_equal(cdf(rel, 0.5, coordinate = 1L), cdf(result, 0.5, coordinate = 1L), tolerance = 1e-8)
})

test_that("sticky Boomerang round-trip preserves segment masks and moments", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  result <- pdmp_sample(
    function(x) x - c(1, -0.5), d = 2,
    flow = "AdaptiveBoomerang", algorithm = "GridThinningStrategy",
    flow_mean = c(1, -0.5), flow_cov = diag(c(2, 0.5)),
    x0 = c(0, 0.4), T = 8, t_warmup = 2,
    sticky = TRUE, can_stick = c(TRUE, TRUE),
    model_prior = bernoulli(0.8),
    slab_prior = dense_gaussian_slab(c(0, 0), diag(2)),
    seed = 94012L, show_progress = FALSE, materialize = FALSE)
  expected_mean <- mean(result)
  expected_pip <- inclusion_probs(result)

  saved <- materialize(result)
  masks <- saved$skeleton[[1L]]$free_masks
  expect_type(masks, "logical")
  expect_equal(dim(masks), dim(saved$skeleton[[1L]]$positions))
  expect_true(any(!masks))
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(saved, tmp)
  restored <- readRDS(tmp)

  expect_equal(mean(restored), expected_mean, tolerance = 1e-10)
  expect_equal(inclusion_probs(restored), expected_pip, tolerance = 1e-10)
  expect_identical(restored$skeleton[[1L]]$free_masks, masks)
})

test_that("sticky ZigZag exports sparse mask transitions and live endpoint", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  result <- pdmp_sample(
    function(x) x, d = 2, flow = "ZigZag",
    algorithm = "GridThinningStrategy",
    x0 = c(0, 0.5), T = 20,
    sticky = TRUE, can_stick = c(TRUE, TRUE),
    model_prior = bernoulli(0.7),
    slab_prior = dense_gaussian_slab(c(0, 0), diag(2)),
    seed = 94013L, show_progress = FALSE, materialize = FALSE)
  saved <- materialize(result)
  sk <- saved$skeleton[[1L]]

  expect_true(isTRUE(sk$sparse))
  expect_type(sk$initial_free_mask, "logical")
  expect_length(sk$initial_free_mask, 2L)
  expect_type(sk$event_free_values, "logical")
  expect_length(sk$event_free_values, length(sk$event_times))
  endpoint <- saved$terminal_states[[1L]]
  expect_length(endpoint$position, 2L)
  expect_length(endpoint$physical_velocity, 2L)
  expect_length(endpoint$stored_frozen_velocity, 2L)
  expect_identical(as.logical(endpoint$free_mask),
    !(endpoint$position == 0 & endpoint$physical_velocity == 0))
})

test_that("typed streaming round-trip stays manifest-only and exact", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()
  directory <- tempfile("typed-stream-")
  dir.create(directory)
  result <- pdmp_sample(function(x) x, d = 3,
    flow = "PreconditionedBPS", algorithm = "GridThinningStrategy",
    T = 5, t_warmup = 1, seed = 94014L, show_progress = FALSE,
    materialize = FALSE, trace_storage = "streaming",
    streaming_trace_directory = directory,
    streaming_trace_buffer_events = 2L)
  expected_mean <- mean(result)
  expected_var <- var(result)
  saved <- materialize(result)
  expect_true(isTRUE(saved$skeleton[[1L]]$streaming))
  expect_null(saved$skeleton[[1L]]$positions)
  expect_true(file.exists(file.path(saved$skeleton[[1L]]$directory,
    "chunk_00000001.bin")))
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(saved, tmp)
  restored <- readRDS(tmp)
  expect_equal(mean(restored), expected_mean, tolerance = 1e-10)
  expect_equal(var(restored), expected_var, tolerance = 1e-10)
})
