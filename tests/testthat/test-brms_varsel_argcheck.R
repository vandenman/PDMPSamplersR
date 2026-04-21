# Unit tests for brms variable selection argument validation and mapping

# ──────────────────────────────────────────────────────────────────────────────
# map_can_stick
# ──────────────────────────────────────────────────────────────────────────────

test_that("map_can_stick marks non-intercept b coefficients as stickable", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "sigma")
  result <- PDMPSamplersR:::map_can_stick(unc_names)
  expect_equal(result, c(FALSE, TRUE, TRUE, FALSE))
})

test_that("map_can_stick honors metadata-derived supported coefficients", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.f1", "sigma")
  result <- PDMPSamplersR:::map_can_stick(
    unc_names,
    supported_coef_names = c("b.x1", "b.x2")
  )
  expect_equal(result, c(FALSE, TRUE, TRUE, FALSE, FALSE))
})

test_that("map_can_stick supports numeric b.<index> unc names with semantic can_stick names", {
  unc_names <- c("b.1", "b.2", "b.3", "Intercept", "sigma")
  result <- PDMPSamplersR:::map_can_stick(
    unc_names,
    supported_coef_names = c("b.x1", "b.x2", "b.x3"),
    user_can_stick = c(x2 = TRUE)
  )
  expect_equal(result, c(FALSE, TRUE, FALSE, FALSE, FALSE))
})

test_that("map_can_stick errors when supported coefficients cannot be aligned", {
  unc_names <- c("b.Intercept", "b.x1", "sigma")
  expect_error(
    PDMPSamplersR:::map_can_stick(
      unc_names,
      supported_coef_names = c("b.x1", "b.x2")
    ),
    "align"
  )
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

test_that("map_can_stick supports named logical vector (short names)", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.x3", "sigma")
  result <- PDMPSamplersR:::map_can_stick(unc_names, user_can_stick = c(x1 = TRUE, x3 = TRUE))
  expect_equal(result, c(FALSE, TRUE, FALSE, TRUE, FALSE))
})

test_that("map_can_stick supports named logical vector (b. prefix)", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "sigma")
  result <- PDMPSamplersR:::map_can_stick(unc_names, user_can_stick = c("b.x2" = TRUE))
  expect_equal(result, c(FALSE, FALSE, TRUE, FALSE))
})

test_that("map_can_stick named vector: unnamed coefficients default to FALSE", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.x3", "sigma")
  result <- PDMPSamplersR:::map_can_stick(unc_names, user_can_stick = c(x1 = TRUE))
  expect_equal(result, c(FALSE, TRUE, FALSE, FALSE, FALSE))
})

test_that("map_can_stick errors on unknown names in named vector", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "sigma")
  expect_error(
    PDMPSamplersR:::map_can_stick(unc_names, user_can_stick = c(x1 = TRUE, z = TRUE)),
    "Unknown"
  )
})

test_that("map_can_stick errors on duplicate names in named vector", {
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "sigma")
  expect_error(
    PDMPSamplersR:::map_can_stick(
      unc_names,
      user_can_stick = c(x1 = TRUE, x1 = FALSE)
    ),
    "duplicate"
  )
})

# ──────────────────────────────────────────────────────────────────────────────
# supported_b_coef_names
# ──────────────────────────────────────────────────────────────────────────────

test_that("supported_b_coef_names accepts simple numeric main effects", {
  df <- data.frame(y = rnorm(8), x1 = rnorm(8), x2 = rnorm(8))
  result <- PDMPSamplersR:::supported_b_coef_names(
    formula = y ~ x1 + x2,
    data = df,
    fe_names = c("x1", "x2")
  )
  expect_equal(result, c("b.x1", "b.x2"))
})

test_that("supported_b_coef_names rejects factor predictors", {
  df <- data.frame(y = rnorm(8), f = factor(rep(c("a", "b"), each = 4)))
  expect_error(
    PDMPSamplersR:::supported_b_coef_names(
      formula = y ~ f,
      data = df,
      fe_names = c("f_b")
    ),
    "numeric"
  )
})

test_that("supported_b_coef_names rejects interactions", {
  df <- data.frame(y = rnorm(8), x1 = rnorm(8), x2 = rnorm(8))
  expect_error(
    PDMPSamplersR:::supported_b_coef_names(
      formula = y ~ x1 * x2,
      data = df,
      fe_names = c("x1", "x2", "x1:x2")
    ),
    "Unsupported terms"
  )
})

test_that("supported_b_coef_names rejects random effects", {
  df <- data.frame(y = rnorm(8), x = rnorm(8), g = factor(rep(1:2, each = 4)))
  expect_error(
    PDMPSamplersR:::supported_b_coef_names(
      formula = y ~ x + (1 | g),
      data = df,
      fe_names = c("x")
    ),
    "Unsupported terms"
  )
})

test_that("supported_b_coef_names rejects one-to-many term expansions", {
  df <- data.frame(y = rnorm(8), x = rnorm(8))
  expect_error(
    PDMPSamplersR:::supported_b_coef_names(
      formula = y ~ x,
      data = df,
      fe_names = c("x", "x_dup")
    ),
    "one fixed-effect coefficient"
  )
})

test_that("supported_b_coef_names supports metadata-only mode", {
  result <- PDMPSamplersR:::supported_b_coef_names(
    fe_names = c("x1", "sx_1", "x1:x2"),
    formula = NULL
  )
  expect_equal(result, c("b.x1", "b.sx_1", "b.x1:x2"))
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

test_that("stickable_coef_names preserves semantic names with numeric b.<index> unc names", {
  unc_names <- c("b.1", "b.2", "Intercept", "sigma")
  result <- PDMPSamplersR:::stickable_coef_names(
    unc_names,
    supported_coef_names = c("b.x1", "b.x2")
  )
  expect_equal(result, c("b.x1", "b.x2"))
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

test_that("derive_parameter_prior uses semantic coefficients for numeric b.<index> unc names", {
  prior <- data.frame(
    prior = c("normal(0, 5)", "normal(0, 10)"),
    class = c("b", "b"),
    coef = c("x1", "x2"),
    group = c("", ""),
    stringsAsFactors = FALSE
  )
  unc_names <- c("b.1", "b.2", "Intercept", "sigma")
  can_stick <- c(TRUE, TRUE, FALSE, FALSE)
  result <- PDMPSamplersR:::derive_parameter_prior(
    prior,
    unc_names,
    can_stick,
    supported_coef_names = c("b.x1", "b.x2")
  )
  expect_equal(result[1], dnorm(0, 0, 5))
  expect_equal(result[2], dnorm(0, 0, 10))
  expect_equal(result[3], 1.0)
  expect_equal(result[4], 1.0)
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
    prior = data.frame(prior = "", class = "", coef = "", group = "", stringsAsFactors = FALSE),
    subsampled = FALSE
  )
  expect_false(result$sticky)
  expect_null(result$can_stick)
})

test_that("validate_brms_sticky allows subsampled", {
  result <- suppressWarnings(PDMPSamplersR:::validate_brms_sticky(
    TRUE, NULL, PDMPSamplersR::bernoulli(0.5), NULL, 3,
    c("b.Intercept", "b.x", "sigma"),
    prior = data.frame(prior = "normal(0, 5)", class = "b", coef = "", group = "", stringsAsFactors = FALSE),
    subsampled = TRUE
  ))
  expect_true(result$sticky)
  expect_length(result$can_stick, 3)
})

test_that("validate_brms_sticky errors without model_prior", {
  expect_error(
    PDMPSamplersR:::validate_brms_sticky(
      TRUE, NULL, NULL, NULL, 3,
      c("b.Intercept", "b.x", "sigma"),
      prior = data.frame(prior = "normal(0, 5)", class = "b", coef = "", group = "", stringsAsFactors = FALSE),
      subsampled = FALSE
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
    TRUE, NULL, PDMPSamplersR::bernoulli(0.5), NULL, 4,
    c("b.Intercept", "b.x1", "b.x2", "sigma"),
    prior = prior, subsampled = FALSE
  )
  expect_true(result$sticky)
  expect_equal(result$can_stick, c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(result$parameter_prior[2], dnorm(0, 0, 5))
  expect_equal(result$parameter_prior[3], dnorm(0, 0, 5))
  expect_equal(result$parameter_prior[1], 1.0)
  expect_equal(result$parameter_prior[4], 1.0)
})

test_that("validate_brms_sticky uses metadata-derived supported coefficients", {
  prior <- data.frame(
    prior = "normal(0, 5)", class = "b", coef = "", group = "",
    stringsAsFactors = FALSE
  )
  result <- PDMPSamplersR:::validate_brms_sticky(
    TRUE, c(x1 = TRUE), PDMPSamplersR::bernoulli(0.5), NULL, 5,
    c("b.Intercept", "b.x1", "b.x2", "b.f1", "sigma"),
    supported_coef_names = c("b.x1", "b.x2"),
    prior = prior, subsampled = FALSE
  )
  expect_equal(result$can_stick, c(FALSE, TRUE, FALSE, FALSE, FALSE))
})

test_that("validate_brms_sticky expands scalar bernoulli prob", {
  prior <- data.frame(
    prior = "normal(0, 5)", class = "b", coef = "", group = "",
    stringsAsFactors = FALSE
  )
  result <- PDMPSamplersR:::validate_brms_sticky(
    TRUE, NULL, PDMPSamplersR::bernoulli(0.5), NULL, 3,
    c("b.Intercept", "b.x", "sigma"),
    prior = prior, subsampled = FALSE
  )
  expect_length(result$model_prior$prob, 3)
})

test_that("validate_brms_sticky accepts user-supplied parameter_prior", {
  prior <- data.frame(
    prior = "", class = "", coef = "", group = "", stringsAsFactors = FALSE
  )
  result <- PDMPSamplersR:::validate_brms_sticky(
    TRUE, NULL, PDMPSamplersR::bernoulli(0.5), c(0.1), 3,
    c("b.Intercept", "b.x", "sigma"),
    prior = prior, subsampled = FALSE
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
      TRUE, NULL, PDMPSamplersR::bernoulli(0.5), c(0.1, 0.2), 3,
      c("b.Intercept", "b.x", "sigma"),
      prior = prior, subsampled = FALSE
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
      TRUE, NULL, PDMPSamplersR::bernoulli(0.5), NULL, 2,
      c("b.Intercept", "sigma"),
      prior = prior, subsampled = FALSE
    ),
    "No supported"
  )
})

# ──────────────────────────────────────────────────────────────────────────────
# build_sticky_inclusion_probs
# ──────────────────────────────────────────────────────────────────────────────

test_that("build_sticky_inclusion_probs keeps names aligned under partial can_stick", {
  incl_raw <- list(chain1 = c(0.0, 0.9, 0.1, 0.8, 0.0))
  can_stick <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.x3", "sigma")

  out <- PDMPSamplersR:::build_sticky_inclusion_probs(incl_raw, can_stick, unc_names)
  expect_equal(out$chain1, c("b.x1" = 0.9, "b.x3" = 0.8))
})

test_that("build_sticky_inclusion_probs uses semantic names with numeric b.<index> unc names", {
  incl_raw <- list(chain1 = c(0.7, 0.2, 0.0, 0.0))
  can_stick <- c(TRUE, TRUE, FALSE, FALSE)
  unc_names <- c("b.1", "b.2", "Intercept", "sigma")

  out <- PDMPSamplersR:::build_sticky_inclusion_probs(
    incl_raw,
    can_stick,
    unc_names,
    supported_coef_names = c("b.x1", "b.x2")
  )
  expect_equal(out$chain1, c("b.x1" = 0.7, "b.x2" = 0.2))
})

test_that("build_sticky_inclusion_probs clamps values to finite probabilities", {
  incl_raw <- list(
    chain1 = c(NaN, 1.3, 0.0, 0.0),
    chain2 = c(-0.2, Inf, 0.0, 0.0)
  )
  can_stick <- c(TRUE, TRUE, FALSE, FALSE)
  unc_names <- c("b.x1", "b.x2", "Intercept", "sigma")

  out <- PDMPSamplersR:::build_sticky_inclusion_probs(incl_raw, can_stick, unc_names)
  expect_equal(out$chain1, c("b.x1" = 0.0, "b.x2" = 1.0))
  expect_equal(out$chain2, c("b.x1" = 0.0, "b.x2" = 0.0))
})

test_that("build_sticky_inclusion_probs errors on length mismatch", {
  incl_raw <- list(chain1 = c(0.0, 0.9, 0.1))
  can_stick <- c(FALSE, TRUE, FALSE, TRUE, FALSE)
  unc_names <- c("b.Intercept", "b.x1", "b.x2", "b.x3", "sigma")

  expect_error(
    PDMPSamplersR:::build_sticky_inclusion_probs(incl_raw, can_stick, unc_names),
    "unexpected length"
  )
})
