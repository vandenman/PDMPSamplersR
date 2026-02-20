# Extracted from test-mvnormal_cmdstan.R:315

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
skip_if_not_installed("BayesFactor")
PDMPSamplersR:::check_for_julia_setup()
paths <- get_stan_paths()
model_path <- file.path(paths$models, "ttest.stan")
skip_if_not(file.exists(model_path), "ttest.stan not found")
set.seed(123)
nx <- 30
ny <- 35
mux <- 0
muy <- 0.5
sigma2 <- 2.1
x1 <- rnorm(nx, mux, sqrt(sigma2))
x2 <- rnorm(ny, muy, sqrt(sigma2))
rscale <- 1.0
analytic_results <- BayesFactor::ttestBF(x = x1, y = x2, rscale = rscale)
analytic_bf <- BayesFactor::extractBF(analytic_results, onlybf = TRUE)
stan_data <- list(nx = nx, ny = ny, x = x1, y = x2, rscale = rscale)
data_path <- tempfile(fileext = ".json")
on.exit(unlink(data_path))
write_stan_json(stan_data, data_path)
result <- pdmp_sample_from_stanmodel(
    model_path, data_path,
    flow = "ZigZag", T = 100000,
    algorithm = "GridThinningStrategy",
    grid_n = 30,
    sticky = TRUE,
    can_stick = c(FALSE, TRUE, FALSE),
    model_prior = bernoulli(0.5),
    parameter_prior = rep(dcauchy(0, 0, rscale), 3),
    show_progress = FALSE
  )
samples <- result$samples
prior_inclusion_prob    <- 0.5
posterior_inclusion_prob <- mean(samples[, 2] != 0)
prior_odds    <- prior_inclusion_prob / (1 - prior_inclusion_prob)
posterior_odds <- posterior_inclusion_prob / (1 - posterior_inclusion_prob)
pdmp_bf <- posterior_odds / prior_odds
expect_equal(sign(log(pdmp_bf)), sign(log(analytic_bf)))
