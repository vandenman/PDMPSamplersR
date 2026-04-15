# Unit tests for brms variable selection argument validation and mapping

# ──────────────────────────────────────────────────────────────────────────────
# map_can_stick
# ──────────────────────────────────────────────────────────────────────────────

test_that("map_can_stick marks non-intercept b coefficients as stickable", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "sigma")
  result <- PDMPSamplersR:::map_can_stick(unc_names)
  expect_equal(result, c(FALSE, TRUE, TRUE, FALSE))
})

test_that("map_can_stick excludes intercept variants", {
  unc_names <- c("b.Intercept", "b.x", "b.z", "sd_group__Intercept")
  result <- PDMPSamplersR:::map_can_stick(unc_names)
  expect_equal(result, c(FALSE, TRUE, TRUE, FALSE))
})

test_that("map_can_stick with no b coefficients returns all FALSE", {
  unc_names <- c("sigma", "shape")
  result <- PDMPSamplersR:::map_can_stick(unc_names)
  expect_equal(result, c(FALSE, FALSE))
})

test_that("map_can_stick handles user_can_stick override", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.x3", "sigma")
  result <- PDMPSamplersR:::map_can_stick(unc_names, user_can_stick = c(TRUE, FALSE, TRUE))
  expect_equal(result, c(FALSE, TRUE, FALSE, TRUE, FALSE))
})

test_that("map_can_stick errors on wrong-length user_can_stick", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "sigma")
  expect_error(
    PDMPSamplersR:::map_can_stick(unc_names, user_can_stick = c(TRUE, FALSE, TRUE)),
    "length"
  )
})

test_that("map_can_stick errors on non-logical user_can_stick", {
  unc_names <- c("b.Intercept", "b.x1", "sigma")
  expect_error(
    PDMPSamplersR:::map_can_stick(unc_names, user_can_stick = c(1, 0)),
    "logical"
  )
})

# ──────────────────────────────────────────────────────────────────────────────
# stickable_coef_names
# ──────────────────────────────────────────────────────────────────────────────

test_that("stickable_coef_names returns correct names", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "sigma")
  result <- PDMPSamplersR:::stickable_coef_names(unc_names)
  expect_equal(result, c("b.x1", "b.x2"))
})

test_that("stickable_coef_names returns empty for no b coefficients", {
  unc_names <- c("sigma", "shape")
  result <- PDMPSamplersR:::stickable_coef_names(unc_names)
  expect_length(result, 0)
})

# ──────────────────────────────────────────────────────────────────────────────
# .parse_and_eval_prior_at_zero
# ──────────────────────────────────────────────────────────────────────────────

test_that("normal(0, s) prior density at zero is correct", {
  result <- PDMPSamplersR:::.parse_and_eval_prior_at_zero("normal(0, 2.5)", "x")
  expect_equal(result, dnorm(0, 0, 2.5))
})

test_that("student_t(df, 0, s) prior density at zero is correct", {
  result <- PDMPSamplersR:::.parse_and_eval_prior_at_zero("student_t(3, 0, 2.5)", "x")
  expect_equal(result, dt(0, df = 3) / 2.5)
})

test_that("non-zero-centered normal prior errors", {
  expect_error(
    PDMPSamplersR:::.parse_and_eval_prior_at_zero("normal(1, 2.5)", "x"),
    "zero-centered"
  )
})

test_that("non-zero-centered student_t prior errors", {
  expect_error(
    PDMPSamplersR:::.parse_and_eval_prior_at_zero("student_t(3, 1, 2.5)", "x"),
    "zero-centered"
  )
})

test_that("unsupported prior family errors", {
  expect_error(
    PDMPSamplersR:::.parse_and_eval_prior_at_zero("cauchy(0, 2.5)", "x"),
    "allowlist"
  )
})

test_that("flat/improper prior errors", {
  expect_error(
    PDMPSamplersR:::.parse_and_eval_prior_at_zero("", "x"),
    "allowlist"
  )
})

# ──────────────────────────────────────────────────────────────────────────────
# derive_parameter_prior
# ──────────────────────────────────────────────────────────────────────────────

test_that("derive_parameter_prior with class-level normal prior", {
  prior <- data.frame(
    prior = c("normal(0, 5)"),
    class = c("b"),
    coef = c(""),
    group = c(""),
    stringsAsFactors = FALSE
  )
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "sigma")
  can_stick <- c(FALSE, TRUE, TRUE, FALSE)
  result <- PDMPSamplersR:::derive_parameter_prior(prior, unc_names, can_stick)
  expect_equal(result[1], 1.0)
  expect_equal(result[2], dnorm(0, 0, 5))
  expect_equal(result[3], dnorm(0, 0, 5))
  expect_equal(result[4], 1.0)
})

test_that("derive_parameter_prior with coefficient-specific priors", {
  prior <- data.frame(
    prior = c("normal(0, 5)", "normal(0, 10)"),
    class = c("b", "b"),
    coef = c("x1", "x2"),
    group = c("", ""),
    stringsAsFactors = FALSE
  )
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "sigma")
  can_stick <- c(FALSE, TRUE, TRUE, FALSE)
  result <- PDMPSamplersR:::derive_parameter_prior(prior, unc_names, can_stick)
  expect_equal(result[2], dnorm(0, 0, 5))
  expect_equal(result[3], dnorm(0, 0, 10))
})

test_that("derive_parameter_prior errors when no prior found", {
  prior <- data.frame(
    prior = character(0), class = character(0),
    coef = character(0), group = character(0),
    stringsAsFactors = FALSE
  )
  unc_names <- c("b.Intercept", "b.x1", "sigma")
  can_stick <- c(FALSE, TRUE, FALSE)
  expect_error(
    PDMPSamplersR:::derive_parameter_prior(prior, unc_names, can_stick),
    "parameter_prior"
  )
})

# ──────────────────────────────────────────────────────────────────────────────
# validate_brms_sticky
# ──────────────────────────────────────────────────────────────────────────────

test_that("validate_brms_sticky returns noop when sticky=FALSE", {
  result <- PDMPSamplersR:::validate_brms_sticky(
    FALSE, NULL, NULL, NULL, 3, c("b.Intercept", "b.x", "sigma"),
    data.frame(prior = "", class = "", coef = "", group = "", stringsAsFactors = FALSE),
    FALSE
  )
  expect_false(result$sticky)
  expect_null(result$can_stick)
})

test_that("validate_brms_sticky allows subsampled", {
  result <- suppressWarnings(PDMPSamplersR:::validate_brms_sticky(
    TRUE, NULL, bernoulli(0.5), NULL, 3,
    c("b.Intercept", "b.x", "sigma"),
    data.frame(prior = "normal(0, 5)", class = "b", coef = "", group = "", stringsAsFactors = FALSE),
    TRUE
  ))
  expect_true(result$sticky)
  expect_length(result$can_stick, 3)
})

test_that("validate_brms_sticky errors without model_prior", {
  expect_error(
    PDMPSamplersR:::validate_brms_sticky(
      TRUE, NULL, NULL, NULL, 3,
      c("b.Intercept", "b.x", "sigma"),
      data.frame(prior = "normal(0, 5)", class = "b", coef = "", group = "", stringsAsFactors = FALSE),
      FALSE
    ),
    "model_prior"
  )
})

test_that("validate_brms_sticky builds correct can_stick and parameter_prior", {
  prior <- data.frame(
    prior = "normal(0, 5)", class = "b", coef = "", group = "",
    stringsAsFactors = FALSE
  )
  result <- PDMPSamplersR:::validate_brms_sticky(
    TRUE, NULL, bernoulli(0.5), NULL, 4,
    c("b.Intercept", "b.x1", "b.x2", "sigma"),
    prior, FALSE
  )
  expect_true(result$sticky)
  expect_equal(result$can_stick, c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(result$parameter_prior[2], dnorm(0, 0, 5))
  expect_equal(result$parameter_prior[3], dnorm(0, 0, 5))
  expect_equal(result$parameter_prior[1], 1.0)
  expect_equal(result$parameter_prior[4], 1.0)
})

test_that("validate_brms_sticky expands scalar bernoulli prob", {
  prior <- data.frame(
    prior = "normal(0, 5)", class = "b", coef = "", group = "",
    stringsAsFactors = FALSE
  )
  result <- PDMPSamplersR:::validate_brms_sticky(
    TRUE, NULL, bernoulli(0.5), NULL, 3,
    c("b.Intercept", "b.x", "sigma"),
    prior, FALSE
  )
  expect_length(result$model_prior$prob, 3)
})

test_that("validate_brms_sticky accepts user-supplied parameter_prior", {
  prior <- data.frame(
    prior = "", class = "", coef = "", group = "", stringsAsFactors = FALSE
  )
  result <- PDMPSamplersR:::validate_brms_sticky(
    TRUE, NULL, bernoulli(0.5), c(0.1), 3,
    c("b.Intercept", "b.x", "sigma"),
    prior, FALSE
  )
  expect_equal(result$parameter_prior[2], 0.1)
  expect_equal(result$parameter_prior[1], 1.0)
  expect_equal(result$parameter_prior[3], 1.0)
})

test_that("validate_brms_sticky errors on wrong-length parameter_prior", {
  prior <- data.frame(
    prior = "", class = "", coef = "", group = "", stringsAsFactors = FALSE
  )
  expect_error(
    PDMPSamplersR:::validate_brms_sticky(
      TRUE, NULL, bernoulli(0.5), c(0.1, 0.2), 3,
      c("b.Intercept", "b.x", "sigma"),
      prior, FALSE
    ),
    "length"
  )
})

test_that("validate_brms_sticky errors when no stickable coefficients", {
  prior <- data.frame(
    prior = "", class = "", coef = "", group = "", stringsAsFactors = FALSE
  )
  expect_error(
    PDMPSamplersR:::validate_brms_sticky(
      TRUE, NULL, bernoulli(0.5), NULL, 2,
      c("b.Intercept", "sigma"),
      prior, FALSE
    ),
    "No supported"
  )
})
