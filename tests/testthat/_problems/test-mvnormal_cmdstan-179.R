# Extracted from test-mvnormal_cmdstan.R:179

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "PDMPSamplersR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
skip_if_no_julia <- function() {
  julia_available <- tryCatch({
    JuliaCall::julia_setup(verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(julia_available, "Julia is not available")
}
get_stan_paths <- function() {
  stan_model_dir <- system.file("stan", "models", package = "PDMPSamplersR")
  stan_data_dir  <- system.file("stan", "data",   package = "PDMPSamplersR")
  if (!nzchar(stan_model_dir)) {
    pkg_root <- testthat::test_path("..", "..")
    stan_model_dir <- file.path(pkg_root, "inst", "stan", "models")
    stan_data_dir  <- file.path(pkg_root, "inst", "stan", "data")
  }
  list(models = stan_model_dir, data = stan_data_dir)
}

# test -------------------------------------------------------------------------
skip_on_cran()
skip_if_no_julia()
PDMPSamplersR:::check_for_julia_setup()
paths <- get_stan_paths()
model_path <- file.path(paths$models, "mvnormal.stan")
skip_if_not(file.exists(model_path), "mvnormal.stan not found")
set.seed(123)
d <- 7
mean_vec   <- rep(0, d)
cov_matrix <- stats::rWishart(1, df = d + 2, Sigma = diag(d))[, , 1]
cov_matrix <- diag(diag(cov_matrix))
stan_data  <- list(N = d, mu = mean_vec, sigma = cov_matrix)
data_path <- tempfile(fileext = ".json")
on.exit(unlink(data_path))
write_stan_json(stan_data, data_path)
mp <- betabernoulli(1, 1)
result <- pdmp_sample_from_stanmodel(
    model_path, data_path,
    flow = "ZigZag", T = 10000,
    algorithm = "GridThinningStrategy",
    sticky = TRUE, can_stick = rep(TRUE, d),
    model_prior = mp,
    parameter_prior = dnorm(rep(0, d), mean_vec, sqrt(diag(cov_matrix))),
    show_progress = FALSE
  )
samples <- result$samples
prior_incl_prob <- rep(mp$a / (mp$a + mp$b), d)
est_incl <- colMeans(samples != 0)
