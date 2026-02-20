# Extracted from test-mvnormal_cmdstan.R:247

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
model_path <- file.path(paths$models, "logistic_regression.stan")
skip_if_not(file.exists(model_path), "logistic_regression.stan not found")
set.seed(123)
n <- 200
d <- 10
d_on <- floor(0.2 * d)
idx_on <- sample(1:d, d_on, FALSE)
intercept <- -3
beta <- rep(0, d)
beta[idx_on] <- rnorm(d_on, 0, 2)
X <- matrix(rnorm(n * d), n, d)
X <- scale(X)
linpred <- intercept + X %*% beta
y <- c(runif(n, 0, 1) <= plogis(linpred))
sd_prior  <- 10
stan_data <- list(N = n, D = d, X = X, y = y, sd_prior = sd_prior)
data_path <- tempfile(fileext = ".json")
on.exit(unlink(data_path))
write_stan_json(stan_data, data_path)
result <- pdmp_sample_from_stanmodel(
    model_path, data_path,
    flow = "ZigZag", T = 20000,
    algorithm = "GridThinningStrategy",
    sticky = TRUE,
    can_stick = c(FALSE, rep(TRUE, d)),
    model_prior = bernoulli(0.5),
    parameter_prior = dnorm(rep(0, d + 1), 0, sd_prior),
    show_progress = FALSE
  )
samples  <- result$samples
est_coef <- colMeans(samples)
