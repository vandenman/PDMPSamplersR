library(PDMPSamplersR)

# Deterministic bivariate, three-wave RI-CLPM data. Each trajectory is the sum
# of a person-specific random intercept and a within-person VAR(1) process.
set.seed(6104)
N <- 36L
P <- 2L
T_wave <- 3L
wave_mean <- c(-0.10, 0.15, 0.00, 0.20, 0.10, 0.25)
A <- matrix(c(0.45, -0.20, 0.25, 0.35), P, P, byrow = TRUE)
Sigma_ri <- matrix(c(0.45, 0.16, 0.16, 0.35), P, P)
Sigma_w1 <- matrix(c(0.55, 0.10, 0.10, 0.50), P, P)
Sigma_e <- matrix(c(0.30, 0.06, 0.06, 0.28), P, P)
rmvn <- function(n, Sigma) {
  matrix(rnorm(n * nrow(Sigma)), n, nrow(Sigma)) %*% chol(Sigma)
}
random_intercept <- rmvn(N, Sigma_ri)
within <- array(0, dim = c(N, P, T_wave))
within[, , 1L] <- rmvn(N, Sigma_w1)
for (wave in 2:T_wave) {
  within[, , wave] <- within[, , wave - 1L] %*% t(A) + rmvn(N, Sigma_e)
}
z <- do.call(cbind, lapply(seq_len(T_wave), function(wave) {
  sweep(random_intercept + within[, , wave], 2L,
        wave_mean[((wave - 1L) * P + 1L):(wave * P)], "+")
}))

slab_mean <- c(0, 0)
slab_cov <- matrix(c(0.20, 0.04, 0.04, 0.16), 2L, 2L)
riclpm_data <- list(
  N = N, P = P, T = T_wave,
  z_sum = colSums(z), z_crossprod = crossprod(z), jitter = 1e-8,
  slab_mean = slab_mean, slab_cov = slab_cov,
  mean_prior_sd = 1.5, lag_prior_sd = 0.5,
  log_sd_prior_mean = -0.4, log_sd_prior_sd = 0.6,
  rho_prior_sd = 0.7
)
riclpm_model <- system.file(
  "stan", "riclpm", "riclpm_sufficient.stan", package = "PDMPSamplersR")
slab <- dense_gaussian_slab(
  slab_mean, slab_cov, coef = "cross_lagged")
common <- list(
  path_to_stanmodel = riclpm_model, standata = riclpm_data,
  algorithm = "GridThinningStrategy", T = 2,
  grid_n = 8L, grid_t_max = 0.15,
  sticky = TRUE, can_stick = "cross_lagged",
  model_prior = betabernoulli(1, 2), slab_prior = slab,
  show_progress = FALSE, materialize = FALSE
)

fits <- lapply(c("ZigZag", "BouncyParticle"), function(flow) {
  do.call(pdmp_sample_from_stanmodel,
          c(common, list(flow = flow, seed = if (flow == "ZigZag") 611L else 612L)))
})
names(fits) <- c("ZigZag", "BouncyParticle")
vapply(fits, function(fit) fit$d, integer(1))
