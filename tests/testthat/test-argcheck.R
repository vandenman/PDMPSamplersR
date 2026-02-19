test_that("validate_type checks type correctly", {
  expect_silent(validate_type(1.0, type = "double", n = 1))
  expect_silent(validate_type(1L,  type = "integer", n = 1))
  expect_silent(validate_type("a", type = "character", n = 1))
  expect_silent(validate_type(TRUE, type = "logical", n = 1))

  expect_error(validate_type(1L, type = "double"), "type")
  expect_error(validate_type(1.0, type = "integer"), "type")
  expect_error(validate_type(1.0, type = "character"), "type")
})

test_that("validate_type checks length correctly", {
  expect_silent(validate_type(c(1.0, 2.0), type = "double", n = 2))
  expect_error(validate_type(c(1.0, 2.0), type = "double", n = 3), "length")
  expect_error(validate_type(1.0, type = "double", n = 2), "length")
})

test_that("validate_type checks dimensions correctly", {
  m <- matrix(1.0, 3, 3)
  expect_silent(validate_type(m, type = "double", dims = c(3, 3)))
  expect_error(validate_type(m, type = "double", dims = c(3, 4)), "dimensions")
  expect_error(validate_type(c(1.0, 2.0), type = "double", dims = c(1, 2)), "dimensions")
})

test_that("validate_type checks positivity correctly", {
  expect_silent(validate_type(1.0, type = "double", n = 1, positive = TRUE))
  expect_error(validate_type(-1.0, type = "double", n = 1, positive = TRUE), "positive")
  expect_error(validate_type(0.0, type = "double", n = 1, positive = TRUE), "positive")
})

test_that("validate_type checks NA correctly", {
  expect_error(validate_type(NA_real_, type = "double"), "missing")
  expect_silent(validate_type(NA_real_, type = "double", allow_missing = TRUE))
})

test_that("validate_type checks Inf correctly", {
  expect_error(validate_type(Inf, type = "double"), "infinite")
  expect_silent(validate_type(Inf, type = "double", allow_inf = TRUE))
})

test_that("cast_integer works correctly", {
  expect_identical(cast_integer(1.0, n = 1), 1L)
  expect_identical(cast_integer(5, n = 1), 5L)
  expect_identical(cast_integer(1.5, n = 1), 1.5)
})

test_that("validate_pdmp_params validates basic parameters", {
  expect_error(validate_pdmp_params(0, "ZigZag", "ThinningStrategy", 100), "positive")
  expect_error(validate_pdmp_params(-1, "ZigZag", "ThinningStrategy", 100), "positive")
  expect_error(validate_pdmp_params(5, "ZigZag", "ThinningStrategy", -100), "positive")
  expect_error(validate_pdmp_params(5, "BadFlow", "ThinningStrategy", 100), "should be one of")
  expect_error(validate_pdmp_params(5, "ZigZag", "BadAlgorithm", 100), "should be one of")
})

test_that("validate_pdmp_params validates t_warmup", {
  expect_error(validate_pdmp_params(5, "ZigZag", "ThinningStrategy", 100, t_warmup = -1), "non-negative")
  expect_error(validate_pdmp_params(5, "ZigZag", "ThinningStrategy", 100, t_warmup = 100), "less than")
  expect_error(validate_pdmp_params(5, "ZigZag", "ThinningStrategy", 100, t_warmup = 150), "less than")
})

test_that("validate_pdmp_params sets defaults correctly", {
  params <- validate_pdmp_params(3, "ZigZag", "ThinningStrategy", 100)
  expect_equal(params$d, 3L)
  expect_equal(params$flow, "ZigZag")
  expect_equal(params$algorithm, "ThinningStrategy")
  expect_equal(params$T, 100)
  expect_equal(params$t0, 0)
  expect_equal(params$t_warmup, 0)
  expect_equal(params$flow_mean, c(0, 0, 0))
  expect_equal(params$flow_cov, diag(3))
  expect_length(params$x0, 3)
  expect_null(params$theta0)
  expect_true(params$show_progress)
  expect_null(params$discretize_dt)
  expect_false(params$sticky)
  expect_equal(params$n_chains, 1L)
  expect_false(params$threaded)
})

test_that("validate_pdmp_params validates flow_cov symmetry", {
  expect_error(
    validate_pdmp_params(2, "ZigZag", "ThinningStrategy", 100,
                         flow_cov = matrix(c(1, 2, 3, 4), 2, 2)),
    "symmetric"
  )
})

test_that("validate_pdmp_params validates theta0 for ZigZag", {
  expect_error(
    validate_pdmp_params(3, "ZigZag", "ThinningStrategy", 100,
                         theta0 = c(1, 0.5, -1)),
    "theta0"
  )
  expect_silent(
    validate_pdmp_params(3, "ZigZag", "ThinningStrategy", 100,
                         theta0 = c(1, -1, 1))
  )
  # BouncyParticle doesn't have the same restriction
  expect_silent(
    validate_pdmp_params(3, "BouncyParticle", "ThinningStrategy", 100,
                         theta0 = c(1, 0.5, -1))
  )
})

test_that("validate_pdmp_params validates sticky parameters", {
  d <- 3
  expect_error(
    validate_pdmp_params(d, "ZigZag", "ThinningStrategy", 100, sticky = TRUE),
    "model_prior"
  )
  expect_error(
    validate_pdmp_params(d, "ZigZag", "ThinningStrategy", 100,
                         sticky = TRUE, model_prior = bernoulli(0.5)),
    "parameter_prior"
  )
  expect_silent(
    validate_pdmp_params(d, "ZigZag", "ThinningStrategy", 100,
                         sticky = TRUE, model_prior = bernoulli(0.5),
                         parameter_prior = rep(1, d),
                         can_stick = rep(TRUE, d))
  )
  expect_silent(
    validate_pdmp_params(d, "ZigZag", "ThinningStrategy", 100,
                         sticky = TRUE, model_prior = betabernoulli(1, 1),
                         parameter_prior = rep(1, d),
                         can_stick = rep(TRUE, d))
  )
})

test_that("validate_pdmp_params validates n_chains and threaded", {
  expect_error(validate_pdmp_params(3, "ZigZag", "ThinningStrategy", 100, n_chains = 0), "positive")
  expect_error(validate_pdmp_params(3, "ZigZag", "ThinningStrategy", 100, n_chains = -1), "positive")
  expect_error(validate_pdmp_params(3, "ZigZag", "ThinningStrategy", 100, threaded = "yes"), "type")
})
