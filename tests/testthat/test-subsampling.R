# R-side logistic regression gradient (no control variate)
# Matches the Julia LogisticRegressionTarget interface:
#   neg gradient = Λ0 (β - μ0)  +  (n/m) X_S' (σ(X_S β) - y_S)
make_logistic_target <- function(X, y, prior_mean, prior_prec) {
  n_obs <- nrow(X)
  d <- ncol(X)
  logistic <- function(z) 1 / (1 + exp(-z))

  grad_sub <- function(beta, indices) {
    prior_part <- prior_prec %*% (beta - prior_mean)
    Xi <- X[indices, , drop = FALSE]
    yi <- y[indices]
    eta <- Xi %*% beta
    pi <- logistic(eta)
    data_part <- t(Xi) %*% (pi - yi)
    scale <- n_obs / length(indices)
    as.vector(prior_part + scale * data_part)
  }

  grad_full <- function(beta) {
    prior_part <- prior_prec %*% (beta - prior_mean)
    eta <- X %*% beta
    p <- logistic(eta)
    data_part <- t(X) %*% (p - y)
    as.vector(prior_part + data_part)
  }

  list(grad_sub = grad_sub, grad_full = grad_full, n_obs = n_obs, d = d)
}

test_that("pdmp_sample_subsampled runs on logistic regression (ZigZag)", {
  skip_on_cran()
  skip_if_no_julia()

  set.seed(2024)
  n <- 200
  d <- 3
  beta_true <- c(0.5, 1.0, -0.8)

  X <- matrix(rnorm(n * (d - 1)), n, d - 1)
  X <- scale(X, center = TRUE, scale = FALSE)
  X <- cbind(1, X)
  eta <- X %*% beta_true
  y <- stats::rbinom(n, 1, 1 / (1 + exp(-eta)))

  prior_mean <- rep(0, d)
  prior_prec <- diag(d)

  target <- make_logistic_target(X, y, prior_mean, prior_prec)

  result <- pdmp_sample_subsampled(
    grad_sub = target$grad_sub,
    n_obs = n, d = d, subsample_size = 20L,
    flow = "ZigZag",
    algorithm = "GridThinningStrategy",
    T = 20000, t_warmup = 0,
    show_progress = FALSE
  )

  expect_s3_class(result, "pdmp_result")
  expect_equal(result$d, d)

  samples <- discretize(result)
  expect_true(is.matrix(samples))
  expect_equal(ncol(samples), d)
  expect_gt(nrow(samples), 50)

  post_mean <- mean(result)
  expect_equal(length(post_mean), d)
  for (j in seq_len(d)) {
    expect_lt(abs(post_mean[j] - beta_true[j]), 0.3)
  }
})

test_that("pdmp_sample_subsampled with full gradient for reflections", {
  skip_on_cran()
  skip_if_no_julia()

  set.seed(42)
  n <- 200
  d <- 3
  beta_true <- c(0.5, 1.0, -0.8)

  X <- matrix(rnorm(n * (d - 1)), n, d - 1)
  X <- scale(X, center = TRUE, scale = FALSE)
  X <- cbind(1, X)
  eta <- X %*% beta_true
  y <- rbinom(n, 1, 1 / (1 + exp(-eta)))

  prior_mean <- rep(0, d)
  prior_prec <- diag(d)

  target <- make_logistic_target(X, y, prior_mean, prior_prec)

  result <- pdmp_sample_subsampled(
    grad_sub = target$grad_sub,
    n_obs = n, d = d, subsample_size = 20L,
    flow = "BouncyParticle",
    algorithm = "GridThinningStrategy",
    T = 20000, t_warmup = 0,
    grad_full = target$grad_full,
    use_full_gradient_for_reflections = TRUE,
    show_progress = FALSE
  )

  expect_s3_class(result, "pdmp_result")

  post_mean <- mean(result)
  for (j in seq_len(d)) {
    expect_lt(abs(post_mean[j] - beta_true[j]), 1.0)
  }
})

test_that("pdmp_sample_subsampled validates inputs", {
  skip_on_cran()

  expect_error(
    pdmp_sample_subsampled(
      grad_sub = "not a function",
      n_obs = 100, d = 2, subsample_size = 10L,
      flow = "ZigZag", algorithm = "ThinningStrategy", T = 100
    ),
    "must be a function"
  )

  expect_error(
    pdmp_sample_subsampled(
      grad_sub = function(x, idx) x,
      n_obs = 10, d = 2, subsample_size = 10L,
      flow = "ZigZag", algorithm = "ThinningStrategy", T = 100
    ),
    "must be less than"
  )

  expect_error(
    pdmp_sample_subsampled(
      grad_sub = function(x, idx) x,
      n_obs = 100, d = 2, subsample_size = 10L,
      flow = "ZigZag", algorithm = "ThinningStrategy", T = 100,
      use_full_gradient_for_reflections = TRUE, grad_full = NULL
    ),
    "grad_full.*must be provided"
  )
})

test_that("pdmp_sample_subsampled matches full gradient pdmp_sample (logistic)", {
  skip_on_cran()
  skip_if_no_julia()

  set.seed(2024)
  n <- 200
  d <- 3
  beta_true <- c(0.5, 1.0, -0.8)

  X <- matrix(rnorm(n * (d - 1)), n, d - 1)
  X <- scale(X, center = TRUE, scale = FALSE)
  X <- cbind(1, X)
  eta <- X %*% beta_true
  y <- rbinom(n, 1, 1 / (1 + exp(-eta)))

  prior_mean <- rep(0, d)
  prior_prec <- diag(d)

  target <- make_logistic_target(X, y, prior_mean, prior_prec)

  full_result <- pdmp_sample(
    f = target$grad_full,
    d = d,
    flow = "AdaptiveBoomerang",
    algorithm = "GridThinningStrategy",
    T = 20000, t_warmup = 0,
    show_progress = FALSE,
    adaptive_scheme = "fullrank"
  )
  full_mean <- mean(full_result)

  sub_result <- pdmp_sample_subsampled(
    grad_sub = target$grad_sub,
    n_obs = n, d = d, subsample_size = as.integer(ceiling(n / 10)),
    flow = "AdaptiveBoomerang",
    algorithm = "GridThinningStrategy",
    T = 20000, t_warmup = 0,
    show_progress = FALSE,
    adaptive_scheme = "fullrank"
  )
  sub_mean <- mean(sub_result)

  for (j in seq_len(d)) {
    expect_lt(abs(full_mean[j] - sub_mean[j]), 0.5)
  }
})
