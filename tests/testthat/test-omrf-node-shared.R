node_shared_omrf_fixture <- function() {
  X <- rbind(
    c(0L, 1L, 0L), c(1L, 0L, 1L), c(1L, 1L, 0L),
    c(0L, 0L, 1L), c(1L, 0L, 0L), c(0L, 1L, 1L)
  )
  P <- ncol(X)
  edges <- do.call(rbind, lapply(seq_len(P - 1L), function(j) {
    cbind(j, seq.int(j + 1L, P))
  }))
  E <- nrow(edges)
  design <- Matrix::sparseMatrix(
    i = rep(seq_len(E), each = 3L),
    j = as.vector(t(cbind(1L, edges[, 1L] + 1L, edges[, 2L] + 1L))),
    x = rep(c(1, 0.5, 0.5), E), dims = c(E, P + 1L)
  )
  data <- list(
    N = nrow(X), P = P, K = 2L, X = X, seen = rep(2L, P),
    prior_only = 0L, prior_threshold_alpha = 2,
    prior_threshold_beta = 2, slab_base_scale = 0.35,
    prior_log_tau_sd = 0.65, prior_log_lambda_sd = 0.75
  )
  list(
    data = data, prior = within(data, prior_only <- 1L), design = design,
    slab = loglinear_gaussian_scale_slab(
      log(0.35), logscale = c("log_tau", "log_lambda"),
      logscale_design = design, coef = "interactions_0")
  )
}

test_that("node-shared OMRF gradient closure retains scale-state priors", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()
  skip_if_not_installed("Matrix")
  fixture <- node_shared_omrf_fixture()
  model <- system.file(
    "stan", "omrf", "omrf_node_shared.stan",
    package = "PDMPSamplersR")
  if (!nzchar(model)) {
    model <- testthat::test_path(
      "..", "..", "inst", "stan", "omrf", "omrf_node_shared.stan")
  }
  compiled <- compile_pdmp_stan_model(model)
  envelope <- omrf_residual_envelope(
    fixture$data$X, fixture$data$seen, "thresholds_0", "interactions_0")
  anchor <- rep(0, 10L)
  subsampling <- stan_subsampling(
    fixture$data$N, 2L, fixture$prior, envelope, anchor = anchor)
  position <- seq(-0.18, 0.27, length.out = 10L)
  subsets <- combn(fixture$data$N, 2L, simplify = FALSE)
  diagnostics <- lapply(subsets, function(subset) {
    stan_subsampling_diagnostics(
      compiled, fixture$data, subsampling, position, subset,
      velocity = seq(0.3, 1.2, length.out = 10L), flow = "ZigZag")
  })
  mean_subsampling <- Reduce(`+`, lapply(diagnostics, `[[`, "subsampled_gradient")) /
    length(diagnostics)
  expect_equal(mean_subsampling, diagnostics[[1L]]$full_gradient, tolerance = 2e-9)
  scale_indices <- 7:10
  expect_gt(max(abs(diagnostics[[1L]]$prior_gradient[scale_indices])), 0.05)
  expect_equal(
    diagnostics[[1L]]$selected_likelihood_gradient[scale_indices],
    rep(0, length(scale_indices)), tolerance = 1e-10)
  expect_identical(fixture$slab$logscale_design$storage, "sparse_rows")
  expect_equal(length(fixture$slab$logscale_design$x), 9L)
})

test_that("public node-shared OMRF samples full and subsampling with ZigZag and BPS", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()
  skip_if_not_installed("Matrix")
  fixture <- node_shared_omrf_fixture()
  model <- system.file(
    "stan", "omrf", "omrf_node_shared.stan",
    package = "PDMPSamplersR")
  if (!nzchar(model)) {
    model <- testthat::test_path(
      "..", "..", "inst", "stan", "omrf", "omrf_node_shared.stan")
  }
  compiled <- compile_pdmp_stan_model(model)
  envelope <- omrf_residual_envelope(
    fixture$data$X, fixture$data$seen, "thresholds_0", "interactions_0")
  subsampling <- stan_subsampling(
    fixture$data$N, 2L, fixture$prior, envelope, anchor = rep(0, 10L))
  common <- list(
    path_to_stanmodel = compiled, standata = fixture$data,
    algorithm = "GridThinningStrategy", T = 1.5,
    grid_n = 6L, grid_t_max = 0.1,
    x0 = c(rep(0, 3L), rep(0.03, 3L), rep(0, 4L)),
    theta0 = c(rep(1, 3L), rep(-1, 3L), rep(1, 4L)),
    sticky = TRUE, can_stick = "interactions_0",
    model_prior = betabernoulli(1, 2), slab_prior = fixture$slab,
    show_progress = FALSE, materialize = FALSE
  )

  for (flow in c("ZigZag", "BouncyParticle")) {
    full_fit <- do.call(pdmp_sample_from_stanmodel, c(
      common, list(flow = flow, seed = if (flow == "ZigZag") 1801L else 1802L)))
    subsampling_fit <- do.call(pdmp_sample_from_stanmodel, c(
      common, list(flow = flow, subsampling = subsampling,
                   seed = if (flow == "ZigZag") 1811L else 1812L)))
    expect_s3_class(full_fit, "pdmp_result")
    expect_s3_class(subsampling_fit, "pdmp_result")
    expect_gt(full_fit$stats$sticky_events[[1L]], 0)
    expect_gt(subsampling_fit$stats$sticky_events[[1L]], 0)
    counters <- attr(subsampling_fit, "subsampling_context_counters")[[1L]]
    if (is.environment(counters)) counters <- as.list(counters)
    expect_equal(
      counters$persons_evaluated,
      subsampling$subsample_size * counters$selected_gradient_calls)
    expect_equal(
      counters$selected_gradient_calls,
      subsampling_fit$stats$residual_oracle_calls[[1L]])
    expect_lte(
      subsampling_fit$stats$subsampling_final_reflections[[1L]],
      subsampling_fit$stats$subsampling_subset_evaluations[[1L]])
    expect_equal(counters$sampling_full_gradient_calls, 0L)
    expect_equal(counters$sampling_model_constructions, 0L)
    expect_equal(counters$sampling_data_constructions, 0L)
  }
})
