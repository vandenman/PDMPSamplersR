test_that("bernoulli creates correct structure", {
  b <- bernoulli(0.5)
  expect_true(is.bernoulli(b))
  expect_false(is.betabernoulli(b))
  expect_equal(b$prob, 0.5)
})

test_that("bernoulli validates inputs", {
  expect_error(bernoulli(-0.1), "between")
  expect_error(bernoulli(1.1), "between")
  expect_error(bernoulli("a"))
})

test_that("bernoulli edge cases", {
  expect_silent(bernoulli(0))
  expect_silent(bernoulli(1))
})

test_that("betabernoulli creates correct structure", {
  bb <- betabernoulli(1, 1)
  expect_true(is.betabernoulli(bb))
  expect_false(is.bernoulli(bb))
  expect_equal(bb$a, 1)
  expect_equal(bb$b, 1)
})

test_that("betabernoulli validates inputs", {
  expect_error(betabernoulli(-1, 1), "positive")
  expect_error(betabernoulli(1, -1), "positive")
})

test_that("is.bernoulli and is.betabernoulli work on arbitrary objects", {
  expect_false(is.bernoulli(list(type = "other")))
  expect_false(is.betabernoulli(list(type = "other")))
  expect_false(is.bernoulli(42))
  expect_false(is.betabernoulli(42))
})

test_that("exchangeable_model_size_prior normalizes positive weights", {
  prior <- exchangeable_model_size_prior(c(1, 2, 1))
  expect_true(PDMPSamplersR:::is.exchangeable_model_size_prior(prior))
  expect_equal(sum(prior$omega), 1)
  expect_equal(prior$omega, c(0.25, 0.5, 0.25))
})

test_that("exchangeable_model_size_prior validates inputs", {
  expect_error(exchangeable_model_size_prior(1), "0:p")
  expect_error(exchangeable_model_size_prior(c(1, 0)), "positive")
  expect_error(exchangeable_model_size_prior(c(1, NA_real_)), "missing")
})
