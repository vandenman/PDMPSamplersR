riclpm_test_data <- function() {
  N <- 10L
  P <- 2L
  T_wave <- 3L
  z <- outer(seq_len(N), seq_len(P * T_wave), function(i, j) {
    sin(0.31 * i + 0.17 * j) + 0.08 * i * j / (N * P * T_wave)
  })
  list(
    N = N, P = P, T = T_wave,
    z_sum = colSums(z), z_crossprod = crossprod(z), jitter = 1e-8,
    slab_mean = c(0.05, -0.03),
    slab_cov = matrix(c(0.20, 0.04, 0.04, 0.16), 2L, 2L),
    mean_prior_sd = 1.5, lag_prior_sd = 0.5,
    log_sd_prior_mean = -0.4, log_sd_prior_sd = 0.6,
    rho_prior_sd = 0.7
  )
}

test_that("bundled RI-CLPM exposes the complete cross-lagged Gaussian slab", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()
  model <- system.file(
    "stan", "riclpm", "riclpm_sufficient.stan", package = "PDMPSamplersR")
  if (!nzchar(model)) {
    model <- testthat::test_path(
      "..", "..", "inst", "stan", "riclpm", "riclpm_sufficient.stan")
  }
  data <- riclpm_test_data()
  mapping <- stan_parameter_mapping(model, data, "cross_lagged")
  expect_equal(nrow(mapping), 2L)
  expect_match(mapping$name, "cross_lagged")

  slab <- dense_gaussian_slab(
    data$slab_mean, data$slab_cov, coef = "cross_lagged")
  for (flow in c("ZigZag", "BouncyParticle", "AdaptiveBoomerang")) {
    adaptive <- identical(flow, "AdaptiveBoomerang")
    fit <- pdmp_sample_from_stanmodel(
      model, data, flow = flow, algorithm = "GridThinningStrategy",
      T = if (adaptive) 10 else 0.2,
      t_warmup = if (adaptive) 2 else 0,
      grid_n = 4L, grid_t_max = if (adaptive) 0.1 else 0.05,
      sticky = TRUE, can_stick = "cross_lagged",
      model_prior = betabernoulli(1, 2), slab_prior = slab,
      show_progress = FALSE, materialize = FALSE,
      seed = match(flow, c("ZigZag", "BouncyParticle",
        "AdaptiveBoomerang")) + 1700L
    )
    expect_s3_class(fit, "pdmp_result")
    expect_equal(fit$d, 19L)
    expect_true(all(is.finite(fit$stats$elapsed_time)))
    expect_true(all(fit$stats$main_events >= 0))
    inclusion <- inclusion_probs(fit)
    expect_length(inclusion, fit$d)
    expect_true(all(is.finite(inclusion)))
    if (adaptive) {
      expect_true(all(fit$stats$warmup_events > 0))
      expect_true(all(fit$stats$main_events > 0))
      expect_true(all(is.finite(fit$stats$final_lambda_ref)))
      expect_true(all(abs(fit$stats$final_lambda_ref - 0.1) > 1e-8))
      expect_true(all(fit$stats$sticky_freezes > 0))
      expect_true(all(fit$stats$sticky_unfreezes > 0))
    }
  }
})
