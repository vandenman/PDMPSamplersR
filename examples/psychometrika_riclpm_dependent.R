library(PDMPSamplersR)

# The RI-CLPM Stan model evaluates its Gaussian likelihood from sufficient
# statistics and contains the complete exchangeable Gaussian slab below.
mapping <- stan_parameter_mapping(
  riclpm_model, riclpm_data, parameters = "cross_lagged")
p <- nrow(mapping)
tau <- 0.35
correlation <- 0.25
slab <- exchangeable_gaussian_slab(
  mean = 0,
  u = tau^2 * (1 - correlation),
  v = tau^2 * correlation,
  coef = "cross_lagged")

fit <- pdmp_sample_from_stanmodel(
  riclpm_model, riclpm_data,
  flow = "BouncyParticle", algorithm = "GridThinningStrategy",
  sticky = TRUE, can_stick = "cross_lagged",
  model_prior = betabernoulli(1, 4), slab_prior = slab)
