test_that("marked Stan specification validates certified envelopes", {
  envelope <- stan_residual_envelope(matrix(c(1, 2, 3), nrow = 1))
  marked <- stan_marked_subsampling(
    n_observations = 3,
    subsample_size = 1,
    prior_standata = list(prior_only = 1L),
    residual_envelope = envelope
  )
  expect_s3_class(marked, "stan_marked_subsampling")
  expect_equal(marked$subsample_size, 1L)
  expect_error(
    stan_marked_subsampling(3, 3, list(), envelope),
    "1 <= m < N"
  )
  expect_error(
    stan_marked_subsampling(3, 1, list(), list(weights = matrix(1, 1, 3))),
    "certified"
  )
})

test_that("OMRF provider uses persons and finite categorical curvature weights", {
  X <- rbind(c(0, 1, 2), c(2, 0, 1), c(1, 2, 0))
  envelope <- omrf_residual_envelope(
    X, seen = c(3, 3, 3),
    thresholds = "thresholds_0",
    interactions = "interactions_0"
  )
  expect_s3_class(envelope, "omrf_residual_envelope")
  expect_equal(ncol(envelope$weights), nrow(X))
  expect_true(all(is.finite(envelope$weights)))
  expect_true(all(envelope$weights >= 0))
  expect_equal(unname(envelope$edge_order),
               rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L)))
  expect_error(
    omrf_residual_envelope(X, c(2, 3, 3), "a", "b"),
    "encoded"
  )
})

test_that("parameter blocks resolve to ordered unconstrained coordinates", {
  names <- c("thresholds_0.1", "thresholds_0.2",
             "interactions_0.1", "log_tau")
  expect_equal(
    PDMPSamplersR:::.resolve_unconstrained_spec(
      c("thresholds_0", "log_tau"), names),
    c(1L, 2L, 4L)
  )
  expect_error(
    PDMPSamplersR:::.resolve_unconstrained_spec(
      c("thresholds_0", "thresholds_0.1"), names),
    "duplicates"
  )
})

test_that("log-linear Gaussian scale slab validates compact design", {
  slab <- loglinear_gaussian_scale_slab(
    log_base_sd = log(0.5),
    logscale = c("log_tau", "log_lambda"),
    logscale_design = matrix(c(1, 0.5, 1, 0.5), nrow = 2),
    coef = "interactions_0"
  )
  expect_s3_class(slab, "dependent_slab_prior")
  expect_equal(slab$type, "loglinear_gaussian_scale")
  expect_error(
    loglinear_gaussian_scale_slab(0, 1L, matrix(1, 2, 2)),
    "Columns"
  )
})

test_that("bundled subset header exposes thread-local install and clear hooks", {
  path <- file.path(testthat::test_path("..", ".."),
                    "inst", "stan", "pdmp_subsample.hpp")
  code <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(code, "thread_local")
  expect_match(code, "pdmp_set_subsample_indices", fixed = TRUE)
  expect_match(code, "pdmp_clear_subsample_indices", fixed = TRUE)
  expect_match(code, "pdmp_get_subsample_index", fixed = TRUE)
})
