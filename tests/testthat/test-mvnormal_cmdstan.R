# Integration tests based on examples/mvnormal_cmdstan.R
# These tests exercise Stan-based sampling with various configurations.

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

# --- 1. MVNormal with ZigZag (flow_mean / flow_cov) -------------------------

test_that("mvnormal ZigZag with known mean and covariance recovers parameters", {
  skip_on_cran()
  skip_if_no_julia()
  PDMPSamplersR:::check_for_julia_setup()

  paths <- get_stan_paths()
  model_path <- file.path(paths$models, "mvnormal.stan")
  skip_if_not(file.exists(model_path), "mvnormal.stan not found")

  set.seed(123)
  d <- 5
  mean_vec   <- rnorm(d, 0, 2)
  cov_matrix <- stats::rWishart(1, df = d + 2, Sigma = diag(d))[, , 1]
  stan_data  <- list(N = d, mu = mean_vec, sigma = cov_matrix)

  data_path <- tempfile(fileext = ".json")
  on.exit(unlink(data_path))
  write_stan_json(stan_data, data_path)

  result <- pdmp_sample_from_stanmodel(
    model_path, data_path,
    flow = "ZigZag", T = 10000,
    flow_mean = mean_vec, flow_cov = cov_matrix,
    show_progress = FALSE
  )

  samples  <- result$samples
  est_mean <- colMeans(samples)
  est_cov  <- cov(samples)

  truth    <- c(mean_vec, cov_matrix[lower.tri(cov_matrix, diag = TRUE)])
  estimate <- c(est_mean, est_cov[lower.tri(est_cov, diag = TRUE)])

  expect_gt(cor(truth, estimate), 0.85)
})

# --- 2. MVNormal with GridThinningStrategy -----------------------------------

test_that("mvnormal ZigZag GridThinningStrategy recovers parameters", {
  skip_on_cran()
  skip_if_no_julia()
  PDMPSamplersR:::check_for_julia_setup()

  paths <- get_stan_paths()
  model_path <- file.path(paths$models, "mvnormal.stan")
  skip_if_not(file.exists(model_path), "mvnormal.stan not found")

  set.seed(123)
  d <- 5
  mean_vec   <- rnorm(d, 0, 2)
  cov_matrix <- stats::rWishart(1, df = d + 2, Sigma = diag(d))[, , 1]
  stan_data  <- list(N = d, mu = mean_vec, sigma = cov_matrix)

  data_path <- tempfile(fileext = ".json")
  on.exit(unlink(data_path))
  write_stan_json(stan_data, data_path)

  result <- pdmp_sample_from_stanmodel(
    model_path, data_path,
    flow = "ZigZag", algorithm = "GridThinningStrategy", T = 10000,
    show_progress = FALSE
  )

  samples  <- result$samples
  est_mean <- colMeans(samples)
  est_cov  <- cov(samples)

  truth    <- c(mean_vec, cov_matrix[lower.tri(cov_matrix, diag = TRUE)])
  estimate <- c(est_mean, est_cov[lower.tri(est_cov, diag = TRUE)])

  expect_gt(cor(truth, estimate), 0.85)
})

# --- 3. Spike-and-slab with diagonal covariance (Bernoulli prior) -----------

test_that("spike-and-slab with diagonal covariance recovers inclusion probabilities", {
  skip_on_cran()
  skip_if_no_julia()
  PDMPSamplersR:::check_for_julia_setup()

  paths <- get_stan_paths()
  model_path <- file.path(paths$models, "mvnormal.stan")
  skip_if_not(file.exists(model_path), "mvnormal.stan not found")

  set.seed(123)
  d <- 12
  mean_vec   <- rep(0, d)
  cov_matrix <- stats::rWishart(1, df = d + 2, Sigma = diag(d))[, , 1]
  cov_matrix <- diag(diag(cov_matrix))
  stan_data  <- list(N = d, mu = mean_vec, sigma = cov_matrix)

  data_path <- tempfile(fileext = ".json")
  on.exit(unlink(data_path))
  write_stan_json(stan_data, data_path)

  prior_incl_prob <- seq(0.1, 0.9, length.out = d)

  result <- pdmp_sample_from_stanmodel(
    model_path, data_path,
    flow = "ZigZag", T = 10000,
    algorithm = "GridThinningStrategy",
    sticky = TRUE, can_stick = rep(TRUE, d),
    model_prior = bernoulli(prior_incl_prob),
    parameter_prior = dnorm(rep(0, d), mean_vec, sqrt(diag(cov_matrix))),
    show_progress = FALSE
  )

  est_incl <- colMeans(result$samples != 0)

  expect_gt(cor(prior_incl_prob, est_incl), 0.85)
  # all estimated inclusion probabilities should be between 0 and 1
  expect_true(all(est_incl >= 0 & est_incl <= 1))
  # MAE should be reasonable (< 0.2)
  expect_lt(mean(abs(prior_incl_prob - est_incl)), 0.2)
})

# --- 4. Spike-and-slab with beta-binomial prior ------------------------------

test_that("spike-and-slab with beta-binomial prior recovers model probabilities", {
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

  # inclusion probabilities: under BetaBernoulli(1,1), marginal inclusion prob is 0.5
  prior_incl_prob <- rep(mp$a / (mp$a + mp$b), d)
  est_incl <- colMeans(samples != 0)
  # expect_gt(cor(prior_incl_prob, est_incl), 0.85) # sd of prior_incl_probs is 0
  expect_lt(sqrt(mean((prior_incl_prob - est_incl)^2)), 0.5)

  # model probabilities: compute observed vs expected
  beta_binomial_pmf <- function(k, d, a, b, log = FALSE) {
    log_beta_num <- lbeta(a + k, b + d - k)
    log_beta_den <- lbeta(a, b)
    log_prob     <- log_beta_num - log_beta_den
    if (log) log_prob else exp(log_prob)
  }

  model_indicators <- apply(1L * (samples != 0), 1L, paste0, collapse = "")
  model_freq  <- table(model_indicators)
  model_probs <- model_freq / sum(model_freq)
  model_sizes <- nchar(gsub("0", "", names(model_freq), fixed = TRUE))

  true_log_probs <- beta_binomial_pmf(model_sizes, d, mp$a, mp$b, log = TRUE)
  est_log_probs  <- log(unname(c(model_probs)))

  expect_gt(cor(true_log_probs, est_log_probs), 0.85)
  # RMSE on the log scale should be moderate (< 1.5)
  expect_lt(sqrt(mean((true_log_probs - est_log_probs)^2)), 1.5)
})

# --- 5. Logistic regression (sd_prior = 10) ----------------------------------

test_that("logistic regression with spike-and-slab recovers parameters", {
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
  true_coef <- c(intercept, beta)

  # correlation between true and estimated coefficients
  expect_gt(cor(true_coef, est_coef), 0.85)

  # active coefficients should have higher inclusion probability than inactive
  p_incl <- colMeans(samples != 0)
  # intercept (col 1) is always included; check betas only
  active_incl   <- mean(p_incl[1 + idx_on])
  inactive_idx  <- setdiff(2:(d + 1), 1 + idx_on)
  inactive_incl <- mean(p_incl[inactive_idx])
  expect_gt(active_incl, inactive_incl)
})

# --- 6. T-test Bayes factor comparison with BayesFactor ----------------------

test_that("t-test Bayes factor is close to analytic result", {
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

  # log Bayes factors should be on the same side (both > 0 or both < 0)
  expect_equal(sign(log(pdmp_bf)), sign(log(analytic_bf)))
  # log BF ratio should be within a factor of 3 on the log scale
  expect_lt(abs(log(pdmp_bf) - log(analytic_bf)), log(3))
})
