test_that("random_effect_subsampling_diagnostics identifies sparse random-effect blocks", {
  N <- 500L
  Z <- matrix(0, nrow = N, ncol = 20)
  group_idx <- split(seq_len(N), rep(seq_len(20), each = 25))
  for (j in seq_len(20)) {
    Z[group_idx[[j]], j] <- 1
  }

  sdata <- list(N = N, Z_1_1 = Z)
  diag <- PDMPSamplersR:::random_effect_subsampling_diagnostics(sdata, 50L)

  expect_equal(nrow(diag), 1L)
  expect_equal(diag$block[[1]], "Z_1_1")
  expect_equal(diag$min_support[[1]], 25)
  expect_equal(diag$expected_support[[1]], 2.5)
  expect_gt(diag$p_zero[[1]], 0.01)
})

test_that("warn_if_low_random_effect_subsampling_support warns for sparse support", {
  N <- 500L
  Z <- matrix(0, nrow = N, ncol = 20)
  group_idx <- split(seq_len(N), rep(seq_len(20), each = 25))
  for (j in seq_len(20)) {
    Z[group_idx[[j]], j] <- 1
  }

  sdata <- list(N = N, Z_1_1 = Z)

  expect_warning(
    PDMPSamplersR:::warn_if_low_random_effect_subsampling_support(sdata, 50L),
    "Subsampled random-effects gradients may be unstable"
  )
})

test_that("warn_if_low_random_effect_subsampling_support is silent for dense support", {
  sdata <- list(N = 500L, Z_1_1 = matrix(1, nrow = 500, ncol = 3))

  expect_no_warning(
    PDMPSamplersR:::warn_if_low_random_effect_subsampling_support(sdata, 50L)
  )
})