# Tests for brm_stancode() and brm_standata() transparency functions.

skip_if_no_brms <- function() {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("rstan")
}

# ── brm_stancode ─────────────────────────────────────────────────────────────

test_that("brm_stancode without subsampling matches brms::stancode", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(20), x = rnorm(20))
  brms_code <- brms::stancode(y ~ x, data = df)
  pdmp_code <- brm_stancode(y ~ x, data = df)
  expect_equal(pdmp_code, brms_code)
})

test_that("brm_stancode with subsampling returns a list with standard and ext_cpp", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, subsample_size = 10L)

  expect_true(is.list(result))
  expect_named(result, c("standard", "ext_cpp"))
  expect_true(is.character(result$standard))
  expect_true(is.character(result$ext_cpp))
})

test_that("brm_stancode standard variant moves means_X to data block", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  original <- brms::stancode(y ~ x, data = df)
  result <- brm_stancode(y ~ x, data = df, subsample_size = 10L)
  modified <- result$standard

  expect_false(identical(original, modified))

  lines <- strsplit(modified, "\n")[[1]]
  data_start <- grep("^data\\s*\\{", lines)
  td_start <- grep("^transformed data\\s*\\{", lines)

  decl_line <- grep("vector\\[Kc\\] means_X", lines)
  expect_length(decl_line, 1)
  expect_true(decl_line > data_start & decl_line < td_start)

  expect_false(any(grepl("means_X\\[i - 1\\] = mean\\(X", modified)))
})

test_that("brm_stancode with subsampling works for bernoulli", {
  skip_if_no_brms()

  df <- data.frame(y = rbinom(50, 1, 0.5), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = brms::bernoulli(),
                          subsample_size = 10L)
  expect_true(is.list(result))
  expect_named(result, c("standard", "ext_cpp"))
  expect_true(grepl("means_X", result$standard))
  expect_true(grepl("get_subsampled_Y_int", result$ext_cpp))
})

# ── brm_standata ─────────────────────────────────────────────────────────────

test_that("brm_standata without subsampling matches brms::standata", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(20), x = rnorm(20))
  brms_data <- brms::standata(y ~ x, data = df)
  pdmp_data <- brm_standata(y ~ x, data = df)
  expect_equal(pdmp_data, brms_data)
})

test_that("brm_standata with subsampling returns three datasets", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_standata(y ~ x, data = df, subsample_size = 10L)

  expect_true(is.list(result))
  expect_named(result, c("full", "prior", "subsample"))
})

test_that("brm_standata full dataset has correct structure", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_standata(y ~ x, data = df, subsample_size = 10L)

  expect_equal(result$full$N, 50L)
  expect_true("means_X" %in% names(result$full))
  expect_equal(length(result$full$means_X), result$full$Kc)
})

test_that("brm_standata prior dataset has N=1 and prior_only=1", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_standata(y ~ x, data = df, subsample_size = 10L)

  expect_equal(result$prior$N, 1L)
  expect_equal(result$prior$prior_only, 1L)
  expect_true("means_X" %in% names(result$prior))
})

test_that("brm_standata subsample has correct size", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_standata(y ~ x, data = df, subsample_size = 10L)

  expect_equal(result$subsample$N, 10L)
  expect_equal(length(result$subsample$Y), 10L)
  expect_equal(nrow(result$subsample$X), 10L)
  expect_true("means_X" %in% names(result$subsample))
})

test_that("brm_standata with explicit indices uses them", {
  skip_if_no_brms()

  df <- data.frame(y = 1:50 + 0.0, x = rnorm(50))
  idx <- c(1L, 5L, 10L, 15L, 20L)
  result <- brm_standata(y ~ x, data = df, subsample_size = 5L, indices = idx)

  expect_equal(result$subsample$N, 5L)
  expect_equal(as.numeric(result$subsample$Y), c(1, 5, 10, 15, 20))
})

test_that("brm_standata rejects mismatched indices length", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  expect_error(
    brm_standata(y ~ x, data = df, subsample_size = 10L, indices = 1:5),
    "indices"
  )
})

test_that("brm_standata rejects subsample_size >= nrow(data)", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(10), x = rnorm(10))
  expect_error(
    brm_standata(y ~ x, data = df, subsample_size = 10L),
    "subsample_size"
  )
})

test_that("brm_standata rejects random effects with subsampling", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(20), x = rnorm(20), g = rep(1:4, each = 5))
  expect_error(
    brm_standata(y ~ x + (1 | g), data = df, subsample_size = 5L),
    "fixed-effects"
  )
})

test_that("brm_standata means_X matches column means of X", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_standata(y ~ x, data = df, subsample_size = 10L)

  sdata <- brms::standata(y ~ x, data = df)
  expected_means <- colMeans(sdata$X[, -1, drop = FALSE])
  expect_equal(as.numeric(result$full$means_X), as.numeric(expected_means))
  expect_equal(as.numeric(result$prior$means_X), as.numeric(expected_means))
  expect_equal(as.numeric(result$subsample$means_X), as.numeric(expected_means))
})

test_that("brm_standata works for bernoulli family", {
  skip_if_no_brms()

  df <- data.frame(y = rbinom(50, 1, 0.5), x = rnorm(50))
  result <- brm_standata(y ~ x, data = df,
                          family = brms::bernoulli(), subsample_size = 10L)

  expect_named(result, c("full", "prior", "subsample"))
  expect_equal(result$subsample$N, 10L)
  expect_equal(result$prior$N, 1L)
})

test_that("brm_standata works for poisson family", {
  skip_if_no_brms()

  df <- data.frame(y = rpois(50, 3), x = rnorm(50))
  result <- brm_standata(y ~ x, data = df,
                          family = poisson(), subsample_size = 10L)

  expect_named(result, c("full", "prior", "subsample"))
  expect_equal(result$subsample$N, 10L)
  expect_true(is.integer(result$prior$Y))
})

test_that("brm_stancode and brm_standata are consistent", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  code <- brm_stancode(y ~ x, data = df, subsample_size = 10L)
  dlist <- brm_standata(y ~ x, data = df, subsample_size = 10L)

  expect_true(grepl("means_X", code$standard))
  expect_true("means_X" %in% names(dlist$full))
  expect_true("means_X" %in% names(dlist$subsample))
})

test_that("brm_stancode ext_cpp injects external C++ declarations for gaussian", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, subsample_size = 10L)

  expect_true(grepl("get_subsampled_Y_real", result$ext_cpp))
  expect_true(grepl("get_subsampled_Xc", result$ext_cpp))
  expect_false(grepl("get_subsampled_Y_real", result$standard))
})

test_that("brm_stancode ext_cpp injects external C++ declarations for poisson", {
  skip_if_no_brms()

  df <- data.frame(y = rpois(50, 3), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = poisson(),
                          subsample_size = 10L)

  expect_true(grepl("get_subsampled_Y_int", result$ext_cpp))
  expect_true(grepl("get_subsampled_Xc", result$ext_cpp))
})
