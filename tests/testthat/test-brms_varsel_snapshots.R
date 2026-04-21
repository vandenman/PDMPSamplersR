# Snapshot tests for fixed-effect variable selection mapping.
#
# These snapshots detect drift in the mapping from unconstrained
# parameter names to can_stick vectors.

test_that("snapshot: can_stick for gaussian y ~ x1 + x2 + x3", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.x3", "sigma")
  result <- PDMPSamplersR:::map_can_stick(unc_names)
  expect_snapshot({
    data.frame(
      parameter = unc_names,
      can_stick = result
    )
  })
})

test_that("snapshot: can_stick for logistic y ~ x1 + x2", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2")
  result <- PDMPSamplersR:::map_can_stick(unc_names)
  expect_snapshot({
    data.frame(
      parameter = unc_names,
      can_stick = result
    )
  })
})

test_that("snapshot: can_stick for gaussian with scale parameter", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.x3", "b.x4", "b.x5", "sigma")
  result <- PDMPSamplersR:::map_can_stick(unc_names)
  expect_snapshot({
    data.frame(
      parameter = unc_names,
      can_stick = result
    )
  })
})

test_that("snapshot: stickable_coef_names for 10-predictor model", {
  unc_names <- c("b.Intercept", paste0("b.x", 1:10), "sigma")
  result <- PDMPSamplersR:::stickable_coef_names(unc_names)
  expect_snapshot(result)
})

test_that("snapshot: can_stick with user override", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.x3", "sigma")
  result <- PDMPSamplersR:::map_can_stick(unc_names, user_can_stick = c(TRUE, FALSE, TRUE))
  expect_snapshot({
    data.frame(
      parameter = unc_names,
      can_stick = result
    )
  })
})

test_that("snapshot: metadata-first mapping for simple numeric formula", {
  df <- data.frame(y = rnorm(8), x1 = rnorm(8), x2 = rnorm(8), f1 = rnorm(8))
  supported <- PDMPSamplersR:::supported_b_coef_names(
    formula = y ~ x1 + x2,
    data = df,
    fe_names = c("x1", "x2")
  )
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.f1", "sigma")
  result <- PDMPSamplersR:::map_can_stick(
    unc_names,
    supported_coef_names = supported
  )

  expect_snapshot({
    data.frame(
      parameter = unc_names,
      can_stick = result
    )
  })
})

test_that("snapshot: metadata-first mapping rejects factor formula", {
  df <- data.frame(y = rnorm(8), f = factor(rep(c("a", "b"), each = 4)))
  expect_snapshot(error = TRUE, {
    PDMPSamplersR:::supported_b_coef_names(
      formula = y ~ f,
      data = df,
      fe_names = c("f_b")
    )
  })
})
