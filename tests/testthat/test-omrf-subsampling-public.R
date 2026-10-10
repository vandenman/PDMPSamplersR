test_that("public independent-slab OMRF supports subsampling AdaptiveBoomerang", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  stan <- system.file(
    "stan", "omrf", "omrf.stan", package = "PDMPSamplersR")
  X <- rbind(
    c(0L, 1L, 0L), c(1L, 0L, 1L), c(1L, 1L, 0L),
    c(0L, 0L, 1L), c(1L, 0L, 0L), c(0L, 1L, 1L)
  )
  full <- list(
    N = 6L, P = 3L, K = 2L, X = X, seen = rep(2L, 3L),
    prior_only = 0L, prior_interaction_sd = 0.45,
    prior_threshold_alpha = 2, prior_threshold_beta = 2)
  prior <- full
  prior$prior_only <- 1L
  compiled <- compile_pdmp_stan_model(stan)
  envelope <- omrf_residual_envelope(
    X, full$seen, "thresholds_0", "interactions_0")
  subsampling <- stan_subsampling(
    6L, 1L, prior, envelope, anchor = rep(0, 6L))
  parameter_prior <- rep(1, 6L)
  parameter_prior[4:6] <- 1 / (sqrt(2 * pi) * full$prior_interaction_sd)
  common <- list(
    path_to_stanmodel = compiled, standata = full,
    flow = "AdaptiveBoomerang", algorithm = "GridThinningStrategy",
    T = 1.5, t_warmup = 0.3, grid_n = 6L, grid_t_max = 0.1,
    x0 = c(rep(0, 3L), rep(0.03, 3L)),
    theta0 = c(rep(1, 3L), rep(-1, 3L)),
    curvature_backend = "finite_difference",
    post_warmup_simplify = TRUE,
    lazy_low_tightness_threshold = 0.2,
    lazy_max_low_tightness_rejections = 2L,
    lazy_max_rejections = 4L,
    sticky = TRUE, can_stick = "interactions_0",
    model_prior = bernoulli(0.35), parameter_prior = parameter_prior,
    show_progress = FALSE, materialize = FALSE
  )
  full_fit <- do.call(
    pdmp_sample_from_stanmodel, c(common, list(seed = 1901L)))
  subsampling_fit <- do.call(
    pdmp_sample_from_stanmodel,
    c(common, list(subsampling = subsampling, seed = 1902L)))
  expect_s3_class(full_fit, "pdmp_result")
  expect_s3_class(subsampling_fit, "pdmp_result")
  expect_gt(full_fit$stats$sticky_events[[1L]], 0)
  expect_gt(subsampling_fit$stats$sticky_events[[1L]], 0)
  counters <- attr(subsampling_fit, "subsampling_context_counters")[[1L]]
  if (is.environment(counters)) counters <- as.list(counters)
  expect_identical(counters$omrf_factorization, "person")
  expect_equal(
    counters$persons_evaluated,
    subsampling$subsample_size * counters$analytic_residual_calls)
  expect_equal(counters$selected_gradient_calls, 0L)
  expect_equal(counters$anchor_cache_gradient_calls, 0L)
  expect_equal(counters$sampling_model_constructions, 0L)
  expect_equal(counters$sampling_data_constructions, 0L)
  expect_equal(counters$sampling_full_gradient_calls, 0L)
  expect_gt(counters$proposal_subset_bound_sum, 0)
  expect_gt(counters$proposal_actual_rate_sum, 0)
  expect_equal(
    counters$proposal_envelope_acceptance_ratio,
    counters$proposal_actual_rate_sum / counters$proposal_subset_bound_sum,
    tolerance = 1e-14)
  expect_equal(
    counters$proposal_certificate_efficiency *
      counters$proposal_intrinsic_acceptance_ratio,
    counters$proposal_envelope_acceptance_ratio,
    tolerance = 1e-12)
})
