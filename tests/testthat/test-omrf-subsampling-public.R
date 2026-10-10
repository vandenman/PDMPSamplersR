test_that("public OMRF diagnostics verify explicit persons and single N/m scaling", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  stan <- system.file(
    "stan", "omrf", "omrf.stan", package = "PDMPSamplersR")
  X <- rbind(c(0L, 1L, 2L), c(2L, 0L, 1L), c(1L, 2L, 0L))
  full <- list(
    N = 3L, P = 3L, K = 3L, X = X, seen = rep(3L, 3L),
    prior_only = 0L, prior_interaction_sd = 0.5,
    prior_threshold_alpha = 2, prior_threshold_beta = 2
  )
  prior <- full
  prior$prior_only <- 1L
  envelope <- omrf_residual_envelope(
    X, full$seen, "thresholds_0", "interactions_0"
  )
  anchor <- rep(0, 9L)
  subsampling <- stan_subsampling(3L, 1L, prior, envelope, anchor = anchor)
  expect_identical(
    subsampling$residual_envelope$analytic_prior$type, "gaussian")
  node_prior <- list(
    N = nrow(X), P = ncol(X), K = 2L, X = X, seen = full$seen,
    prior_only = 1L, prior_threshold_alpha = 2,
    prior_threshold_beta = 2, slab_base_scale = 0.7,
    prior_global_scale_sd = 1.2, prior_node_scale_sd = 0.8)
  node_subsampling <- stan_subsampling(
    nrow(X), 1L, node_prior, envelope, anchor = anchor)
  expect_identical(
    node_subsampling$residual_envelope$analytic_prior$type,
    "node_scale_gaussian")
  expect_identical(
    node_subsampling$residual_envelope$analytic_prior$global_scale,
    "global_scale")
  expect_identical(
    node_subsampling$residual_envelope$analytic_prior$node_scale,
    "node_scale")

  compiled <- compile_pdmp_stan_model(stan)
  position <- seq(-0.21, 0.27, length.out = 9L)
  diagnostics <- lapply(seq_len(3L), function(person) {
    stan_subsampling_diagnostics(
      compiled, full, subsampling, position = position, subset = person,
      velocity = rep(1, 9L), flow = "ZigZag")
  })
  # With one person per proposal, N/m = 3 is applied exactly once: the
  # average over persons of the subsampled gradient is the full gradient.
  mean_subsampled <- Reduce(`+`, lapply(diagnostics, `[[`,
    "subsampled_gradient")) / length(diagnostics)
  expect_equal(mean_subsampled, diagnostics[[1L]]$full_gradient,
               tolerance = 1e-9)
  # The residual also subtracts the anchor Hessian term, so it is not the
  # plain first-order difference; its rate is certified per person.
  for (diagnostic in diagnostics) {
    expect_lte(diagnostic$residual_rate,
               diagnostic$envelope_rate * (1 + 1e-10) + 1e-10)
  }

  anchor_diagnostic <- stan_subsampling_diagnostics(
    compiled, full, subsampling, position = anchor, subset = 1L,
    velocity = rep(1, 9L), flow = "BouncyParticle")
  expect_equal(
    anchor_diagnostic$subsampled_gradient, anchor_diagnostic$full_gradient,
    tolerance = 1e-10
  )
  expect_equal(anchor_diagnostic$anchor_closure_error, 0, tolerance = 1e-12)
})

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
