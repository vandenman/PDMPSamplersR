library(PDMPSamplersR)

# Small deterministic ordinal data; persons are the marked observations.
set.seed(2718)
N <- 30L
P <- 4L
K <- 3L
X <- matrix(sample.int(K, N * P, replace = TRUE) - 1L, N, P)
stan_data <- list(
  N = N, P = P, K = K, X = X, seen = rep(K, P), prior_only = 0L,
  prior_interaction_sd = 0.45,
  prior_threshold_alpha = 2,
  prior_threshold_beta = 2
)
prior_data <- stan_data
prior_data$prior_only <- 1L
stan_file <- system.file("stan", "omrf", "omrf_marked.stan",
                         package = "PDMPSamplersR")
compiled <- compile_pdmp_stan_model(stan_file)

slab <- independent_slab_density(
  1 / (sqrt(2 * pi) * stan_data$prior_interaction_sd),
  coef = "interactions_0")
common <- list(
  path_to_stanmodel = compiled, standata = stan_data,
  flow = "ZigZag", algorithm = "GridThinningStrategy",
  T = 50, grid_n = 12L, grid_t_max = 0.5,
  sticky = TRUE, can_stick = "interactions_0",
  model_prior = betabernoulli(1, 4), slab_prior = slab,
  show_progress = FALSE
)

full_fit <- do.call(
  pdmp_sample_from_stanmodel,
  c(common, list(seed = 812))
)

envelope <- omrf_residual_envelope(
  stan_data$X, stan_data$seen,
  thresholds = "thresholds_0",
  interactions = "interactions_0")
marked <- stan_marked_subsampling(stan_data$N, 6L, prior_data, envelope)

marked_fit <- do.call(
  pdmp_sample_from_stanmodel,
  c(common, list(marked_subsampling = marked, seed = 813))
)

rbind(full = mean(full_fit), marked = mean(marked_fit))
