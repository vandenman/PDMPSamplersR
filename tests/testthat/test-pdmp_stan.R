# Integration tests for Stan model sampling.
# These tests are skipped if Julia or the Stan models are not available.

# Pure R-side validation tests (don't call check_for_julia_setup)

test_that("pdmp_sample_from_stanmodel validates input types", {
  expect_error(pdmp_sample_from_stanmodel(123, "data.json"), "type")
  expect_error(pdmp_sample_from_stanmodel("model.stan", 123), "type")
})

test_that("pdmp_sample_from_stanmodel accepts standata as list", {
  model_file <- tempfile(fileext = ".stan")
  file.create(model_file)
  on.exit(unlink(model_file), add = TRUE)

  captured <- new.env(parent = emptyenv())
  captured$data <- NULL
  captured$file <- NULL

  testthat::local_mocked_bindings(
    write_stan_json = function(data, file, always_decimal = FALSE) {
      captured$data <- data
      captured$file <- file
      jsonlite::write_json(list(N = 1), path = file, auto_unbox = TRUE)
    },
    check_for_julia_setup = function() {
      stop("SENTINEL_CHECK_SETUP", call. = FALSE)
    },
    .package = "PDMPSamplersR"
  )

  expect_error(
    pdmp_sample_from_stanmodel(model_file, list(N = 1)),
    "SENTINEL_CHECK_SETUP"
  )

  expect_equal(captured$data, list(N = 1))
  expect_true(is.character(captured$file))
  expect_match(captured$file, "\\.json$")
  expect_false(file.exists(captured$file))
})

test_that("pdmp_sample_from_stanmodel validates file existence", {
  expect_error(pdmp_sample_from_stanmodel("nonexistent.stan", "data.json"), "not found")
  expect_error(pdmp_sample_from_stanmodel("nonexistent.so",   "data.json"), "not found")
})

test_that("pdmp_sample_from_stanmodel gates dependent slabs pending target composition", {
  model_file <- tempfile(fileext = ".stan")
  data_file <- tempfile(fileext = ".json")
  file.create(model_file)
  jsonlite::write_json(list(N = 1), data_file, auto_unbox = TRUE)
  on.exit({unlink(model_file); unlink(data_file)}, add = TRUE)

  expect_error(
    pdmp_sample_from_stanmodel(
      model_file,
      data_file,
      sticky = TRUE,
      algorithm = "GridThinningStrategy",
      can_stick = TRUE,
      model_prior = bernoulli(0.5),
      slab_prior = dense_gaussian_slab(0, matrix(1, 1, 1), coef = 1)
    ),
    "target composition|temporarily gated"
  )
})

test_that("dependent slab Stan validation happens before Julia setup", {
  model_file <- tempfile(fileext = ".stan")
  data_file <- tempfile(fileext = ".json")
  prior_data_file <- tempfile(fileext = ".json")
  file.create(model_file)
  jsonlite::write_json(list(N = 1), data_file, auto_unbox = TRUE)
  jsonlite::write_json(list(N = 1), prior_data_file, auto_unbox = TRUE)
  on.exit(unlink(c(model_file, data_file, prior_data_file)), add = TRUE)

  testthat::local_mocked_bindings(
    check_for_julia_setup = function() {
      stop("SENTINEL_CHECK_SETUP", call. = FALSE)
    },
    .package = "PDMPSamplersR"
  )

  expect_error(
    pdmp_sample_from_stanmodel(
      model_file,
      data_file,
      prior_standata = prior_data_file,
      sticky = FALSE,
      model_prior = bernoulli(0.5),
      slab_prior = dense_gaussian_slab(0, matrix(1, 1, 1), coef = 1)
    ),
    "sticky"
  )

  expect_error(
    pdmp_sample_from_stanmodel(
      model_file,
      data_file,
      prior_standata = prior_data_file,
      sticky = TRUE,
      model_prior = bernoulli(0.5),
      slab_prior = dense_gaussian_slab(0, matrix(1, 1, 1), coef = 1)
    ),
    "GridThinningStrategy"
  )
})

test_that("pdmp_sample_from_stanmodel rejects wrong file extensions", {
  # Create temp files with wrong extensions to test extension check
  wrong_ext <- tempfile(fileext = ".txt")
  data_file <- tempfile(fileext = ".json")
  file.create(wrong_ext)
  file.create(data_file)
  on.exit({unlink(wrong_ext); unlink(data_file)})

  expect_error(pdmp_sample_from_stanmodel(wrong_ext, data_file), "path_to_stanmodel")
})

test_that("pdmp_sample_from_stanmodel rejects wrong data extension", {
  model_file <- tempfile(fileext = ".stan")
  wrong_data <- tempfile(fileext = ".csv")
  file.create(model_file)
  file.create(wrong_data)
  on.exit({unlink(model_file); unlink(wrong_data)})

  expect_error(pdmp_sample_from_stanmodel(model_file, wrong_data), "JSON")
})

test_that("pdmp_sample_from_stanmodel gates one-dimensional named independent slab target", {
  skip_on_cran()

  stan_file <- tempfile(fileext = ".stan")
  data_file <- tempfile(fileext = ".json")
  prior_data_file <- tempfile(fileext = ".json")
  on.exit(unlink(c(stan_file, data_file, prior_data_file)), add = TRUE)

  writeLines(c(
    "data {",
    "  int<lower=0,upper=1> prior_only;",
    "}",
    "parameters {",
    "  real beta;",
    "}",
    "model {",
    "  beta ~ normal(0, 1);",
    "  if (!prior_only) beta ~ normal(0.25, 1);",
    "}"
  ), stan_file)
  write_stan_json(list(prior_only = 0L), data_file)
  write_stan_json(list(prior_only = 1L), prior_data_file)

  expect_error(
    pdmp_sample_from_stanmodel(
      stan_file,
      data_file,
      prior_standata = prior_data_file,
      flow = "ZigZag",
      algorithm = "GridThinningStrategy",
      T = 1.0,
      x0 = 0.2,
      theta0 = 1,
      grid_n = 3L,
      grid_t_max = 0.1,
      sticky = TRUE,
      can_stick = TRUE,
      model_prior = bernoulli(0.5),
      slab_prior = independent_slab_density(1, coef = "beta"),
      show_progress = FALSE,
      materialize = FALSE
    ),
    "target composition|temporarily gated"
  )
})

test_that("pdmp_sample_from_stanmodel runs with mvnormal Stan model", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  stan_model_dir <- system.file("stan", "models", package = "PDMPSamplersR")
  stan_data_dir  <- system.file("stan", "data",   package = "PDMPSamplersR")

  # Fall back to package source directory if not installed
  if (!nzchar(stan_model_dir)) {
    pkg_root <- testthat::test_path("..", "..")
    stan_model_dir <- file.path(pkg_root, "inst", "stan", "models")
    stan_data_dir  <- file.path(pkg_root, "inst", "stan", "data")
  }

  model_path <- file.path(stan_model_dir, "mvnormal.stan")
  skip_if_not(file.exists(model_path), "mvnormal.stan not found")

  d <- 3
  mu <- rep(0, d)
  sigma <- diag(d)
  data_path <- tempfile(fileext = ".json")
  on.exit(unlink(data_path))
  write_stan_json(list(N = d, mu = mu, sigma = sigma), data_path)

  result <- pdmp_sample_from_stanmodel(
    model_path, data_path,
    flow = "ZigZag", T = 1000,
    flow_mean = mu, flow_cov = sigma,
    show_progress = FALSE
  )

  expect_s3_class(result, "pdmp_result")
  expect_equal(result$d, d)

  samples <- discretize(result)
  expect_true(is.matrix(samples))
  expect_equal(ncol(samples), d)
  expect_gt(nrow(samples), 10)

  # posterior means should be near zero for standard normal
  posterior_means <- mean(result)
  expect_true(all(abs(posterior_means) < 1.0))
})
