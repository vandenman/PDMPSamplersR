library(PDMPSamplersR)

# `stan_data` contains N, P, K, X, seen, prior_only = 0, and prior scales.
# `prior_data` is identical except for prior_only = 1.
stan_file <- system.file("stan", "omrf", "omrf_marked.stan",
                         package = "PDMPSamplersR")
compiled <- compile_pdmp_stan_model(stan_file)

envelope <- omrf_residual_envelope(
  stan_data$X, stan_data$seen,
  thresholds = "thresholds_0",
  interactions = "interactions_0")
marked <- stan_marked_subsampling(
  stan_data$N, min(50L, stan_data$N - 1L), prior_data, envelope)

fit <- pdmp_sample_from_stanmodel(
  compiled, stan_data, marked_subsampling = marked,
  flow = "ZigZag", algorithm = "GridThinningStrategy",
  sticky = TRUE, can_stick = "interactions_0",
  model_prior = betabernoulli(1, 4),
  slab_prior = independent_slab_density(
    1 / (sqrt(2 * pi) * stan_data$prior_interaction_sd),
    coef = "interactions_0"))
