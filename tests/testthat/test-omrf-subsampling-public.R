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
  legacy_node_sum <- function(x) {
    edges <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
    total <- 0
    for (j in seq_len(3L)) {
      B <- matrix(0, 2L, 9L)
      for (u in seq_len(2L)) {
        B[u, 2L * (j - 1L) + u] <- 1
        incident <- which(edges[, 1L] == j | edges[, 2L] == j)
        neighbours <- ifelse(
          edges[incident, 1L] == j,
          edges[incident, 2L], edges[incident, 1L])
        B[u, 6L + incident] <- u * x[neighbours]
      }
      total <- total + 0.5 * norm(B, type = "2")^2
    }
    total
  }
  legacy_weights <- apply(X, 1L, legacy_node_sum)
  stacked_weight <- function(x) {
    edges <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
    B <- matrix(0, 6L, 9L)
    for (j in seq_len(3L)) {
      for (u in seq_len(2L)) {
        row <- 2L * (j - 1L) + u
        B[row, row] <- 1
        incident <- which(edges[, 1L] == j | edges[, 2L] == j)
        neighbours <- ifelse(edges[incident, 1L] == j,
                             edges[incident, 2L], edges[incident, 1L])
        B[row, 6L + incident] <- u * x[neighbours]
      }
    }
    0.5 * norm(B, type = "2")^2
  }
  expect_equal(drop(envelope$weights), apply(X, 1L, stacked_weight),
               tolerance = 1e-14)
  duplicate_envelope <- omrf_residual_envelope(
    rbind(X, X[2L, ]), full$seen, "thresholds_0", "interactions_0")
  expect_identical(drop(duplicate_envelope$weights[, c(2L, 4L)]),
                   rep(envelope$weights[, 2L], 2L))
  expect_identical(envelope$bound_type, "stacked_person_spectral")
  expect_true(all(drop(envelope$weights) <= legacy_weights * (1 + 1e-12)))
  expect_true(any(drop(envelope$weights) < legacy_weights * (1 - 1e-8)))
  anchor <- rep(0, 9L)
  subsampling <- stan_subsampling(3L, 2L, prior, envelope, anchor = anchor)
  compiled <- compile_pdmp_stan_model(stan)
  position <- seq(-0.21, 0.27, length.out = 9L)
  subset <- c(1L, 3L)
  diagnostic <- stan_subsampling_diagnostics(
    compiled, full, subsampling, position = position, subset = subset,
    velocity = rep(1, 9L), flow = "ZigZag"
  )
  explicit_person_gradient <- function(q, x) {
    thresholds <- matrix(q[1:6], nrow = 3L, byrow = TRUE)
    interactions <- matrix(0, 3L, 3L)
    interactions[upper.tri(interactions)] <- q[7:9]
    interactions <- interactions + t(interactions)
    field <- drop(x %*% interactions)
    probabilities <- lapply(seq_len(3L), function(j) {
      eta <- c(0, thresholds[j, ] + seq_len(2L) * field[j])
      value <- exp(eta - max(eta))
      value / sum(value)
    })
    result <- numeric(9L)
    for (j in seq_len(3L)) {
      for (u in seq_len(2L)) {
        result[2L * (j - 1L) + u] <-
          probabilities[[j]][u + 1L] - as.numeric(x[j] == u)
      }
    }
    edges <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L))
    for (edge in seq_len(3L)) {
      j <- edges[edge, 1L]
      k <- edges[edge, 2L]
      expected_j <- sum(seq_len(2L) * probabilities[[j]][2:3])
      expected_k <- sum(seq_len(2L) * probabilities[[k]][2:3])
      result[6L + edge] <- -(2 * x[j] * x[k] -
        x[k] * expected_j - x[j] * expected_k)
    }
    result
  }
  explicit_selected <- Reduce(
    `+`, lapply(subset, function(n) explicit_person_gradient(position, X[n, ])))
  expect_equal(
    diagnostic$selected_likelihood_gradient, explicit_selected,
    tolerance = 2e-9)
  expect_equal(diagnostic$selected_scale, 3 / 2)
  expect_equal(
    diagnostic$subsampled_gradient,
    diagnostic$deterministic_gradient +
      diagnostic$selected_scale * diagnostic$residual_gradient,
    tolerance = 1e-12)
  expect_gt(max(abs(diagnostic$residual_gradient)), 1e-4)

  diagnostic_velocity <- diagnostic$residual_gradient
  for (flow in c("ZigZag", "BouncyParticle", "AdaptiveBoomerang")) {
    flow_diagnostic <- stan_subsampling_diagnostics(
      compiled, full, subsampling, position = position, subset = subset,
      velocity = diagnostic_velocity, flow = flow)
    expect_gt(flow_diagnostic$residual_rate, 0)
    expect_lte(
      flow_diagnostic$residual_rate,
      flow_diagnostic$envelope_rate * (1 + 1e-10) + 1e-10)
  }

  anchor_diagnostic <- stan_subsampling_diagnostics(
    compiled, full, subsampling, position = anchor, subset = subset,
    velocity = rep(1, 9L), flow = "BouncyParticle")
  expect_equal(
    anchor_diagnostic$subsampled_gradient, anchor_diagnostic$full_gradient,
    tolerance = 1e-10
  )
  expect_equal(anchor_diagnostic$anchor_closure_error, 0, tolerance = 1e-12)
  expect_equal(diagnostic$model_constructions, 3L)
  expect_equal(diagnostic$data_constructions, 3L)
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
    6L, 2L, prior, envelope, anchor = rep(0, 6L))
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
  expect_equal(
    counters$persons_evaluated,
    subsampling$subsample_size * counters$selected_gradient_calls)
  expect_equal(counters$sampling_model_constructions, 0L)
  expect_equal(counters$sampling_data_constructions, 0L)
  expect_equal(counters$sampling_full_gradient_calls, 0L)
})
