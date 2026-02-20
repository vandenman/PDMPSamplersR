# Integration tests for Stan model sampling.
# These tests are skipped if Julia or the Stan models are not available.

skip_if_no_julia <- function() {
  julia_available <- tryCatch({
    JuliaCall::julia_setup(verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(julia_available, "Julia is not available")

  # pdmp_available <- tryCatch({
  #   PDMPSamplersR:::check_for_julia_setup()
  #   # check BridgeStan extension loaded: PDMPModel(String, String) must exist
  #   JuliaCall::julia_eval("hasmethod(PDMPModel, Tuple{String, String})")
  # }, error = function(e) FALSE)
  # testthat::skip_if_not(pdmp_available, "PDMPSamplers.jl BridgeStan extension not available")
}

# Pure R-side validation tests (don't call check_for_julia_setup)

test_that("pdmp_sample_from_stanmodel validates file path types", {
  expect_error(pdmp_sample_from_stanmodel(123, "data.json"), "type")
  expect_error(pdmp_sample_from_stanmodel("model.stan", 123), "type")
})

test_that("pdmp_sample_from_stanmodel validates file existence", {
  expect_error(pdmp_sample_from_stanmodel("nonexistent.stan", "data.json"), "not found")
  expect_error(pdmp_sample_from_stanmodel("nonexistent.so",   "data.json"), "not found")
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

test_that("pdmp_sample_from_stanmodel runs with mvnormal Stan model", {
  skip_on_cran()
  skip_if_no_julia()
  PDMPSamplersR:::check_for_julia_setup()

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
