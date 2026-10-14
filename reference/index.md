# Package index

## Sampling

- [`pdmp_sample()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmp_sample.md)
  : PDMP Sampling
- [`pdmp_sample_from_stanmodel()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmp_sample_from_stanmodel.md)
  : PDMP Sampling from Stan Model
- [`sampler_control()`](https://vandenman.github.io/PDMPSamplersR/reference/sampler_control.md)
  : Additional sampler controls
- [`stan_compile_options()`](https://vandenman.github.io/PDMPSamplersR/reference/stan_compile_options.md)
  : Stan compilation options
- [`support_boundary_control()`](https://vandenman.github.io/PDMPSamplersR/reference/support_boundary_control.md)
  : Support Boundary Control
- [`write_stan_json()`](https://vandenman.github.io/PDMPSamplersR/reference/write_stan_json.md)
  : Write data to a JSON file readable by Stan
- [`pdmp_subsample_hpp_path()`](https://vandenman.github.io/PDMPSamplersR/reference/stan_subsampling.md)
  [`compile_pdmp_stan_model()`](https://vandenman.github.io/PDMPSamplersR/reference/stan_subsampling.md)
  [`stan_residual_envelope()`](https://vandenman.github.io/PDMPSamplersR/reference/stan_subsampling.md)
  [`stan_subsampling()`](https://vandenman.github.io/PDMPSamplersR/reference/stan_subsampling.md)
  [`omrf_residual_envelope()`](https://vandenman.github.io/PDMPSamplersR/reference/stan_subsampling.md)
  [`stan_parameter_mapping()`](https://vandenman.github.io/PDMPSamplersR/reference/stan_subsampling.md)
  : Compile and configure subsampling custom-Stan models
- [`stan_subsampling_diagnostics()`](https://vandenman.github.io/PDMPSamplersR/reference/stan_subsampling_diagnostics.md)
  : Diagnose a subsampling custom-Stan gradient and envelope
- [`position_anchor_bank()`](https://vandenman.github.io/PDMPSamplersR/reference/position_anchor_bank.md)
  : Frozen position anchor bank from a set of states
- [`online_summary_options()`](https://vandenman.github.io/PDMPSamplersR/reference/online_summary_options.md)
  : Online summaries of the retained trajectory

## Results

- [`materialize()`](https://vandenman.github.io/PDMPSamplersR/reference/materialize.md)
  : Materialise a PDMP result for serialisation
- [`read_online_summaries()`](https://vandenman.github.io/PDMPSamplersR/reference/read_online_summaries.md)
  : Read online summaries written by a streaming chain
- [`mean(`*`<pdmp_result>`*`)`](https://vandenman.github.io/PDMPSamplersR/reference/mean.pdmp_result.md)
  : Continuous-time mean of PDMP trace
- [`median(`*`<pdmp_result>`*`)`](https://vandenman.github.io/PDMPSamplersR/reference/median.pdmp_result.md)
  : Continuous-time median of PDMP trace
- [`var(`*`<pdmp_result>`*`)`](https://vandenman.github.io/PDMPSamplersR/reference/var.pdmp_result.md)
  : Continuous-time variance of PDMP trace
- [`var()`](https://vandenman.github.io/PDMPSamplersR/reference/var.md)
  : Variance
- [`sd()`](https://vandenman.github.io/PDMPSamplersR/reference/sd.md) :
  Standard deviation
- [`cov(`*`<pdmp_result>`*`)`](https://vandenman.github.io/PDMPSamplersR/reference/cov.pdmp_result.md)
  : Continuous-time covariance of PDMP trace
- [`cov()`](https://vandenman.github.io/PDMPSamplersR/reference/cov.md)
  : Covariance matrix
- [`cor(`*`<pdmp_result>`*`)`](https://vandenman.github.io/PDMPSamplersR/reference/cor.pdmp_result.md)
  : Continuous-time correlation of PDMP trace
- [`cor()`](https://vandenman.github.io/PDMPSamplersR/reference/cor.md)
  : Correlation matrix
- [`quantile(`*`<pdmp_result>`*`)`](https://vandenman.github.io/PDMPSamplersR/reference/quantile.pdmp_result.md)
  : Continuous-time quantile of PDMP trace
- [`ess()`](https://vandenman.github.io/PDMPSamplersR/reference/ess.md)
  : Effective sample size
- [`inclusion_probs()`](https://vandenman.github.io/PDMPSamplersR/reference/inclusion_probs.md)
  : Marginal inclusion probabilities
- [`cdf()`](https://vandenman.github.io/PDMPSamplersR/reference/cdf.md)
  : Empirical CDF
- [`discretize()`](https://vandenman.github.io/PDMPSamplersR/reference/discretize.md)
  : Discretize a PDMP trace

## Transforms

- [`identity_transform()`](https://vandenman.github.io/PDMPSamplersR/reference/identity_transform.md)
  : Identity transform (no constraint)
- [`lower_transform()`](https://vandenman.github.io/PDMPSamplersR/reference/lower_transform.md)
  : Lower-bounded transform
- [`upper_transform()`](https://vandenman.github.io/PDMPSamplersR/reference/upper_transform.md)
  : Upper-bounded transform
- [`double_transform()`](https://vandenman.github.io/PDMPSamplersR/reference/double_transform.md)
  : Double-bounded transform

## Priors

- [`bernoulli()`](https://vandenman.github.io/PDMPSamplersR/reference/bernoulli.md)
  : Create a Bernoulli model prior
- [`betabernoulli()`](https://vandenman.github.io/PDMPSamplersR/reference/betabernoulli.md)
  : Create a Beta-Bernoulli model prior
- [`dense_gaussian_slab()`](https://vandenman.github.io/PDMPSamplersR/reference/dependent_slab_priors.md)
  [`exchangeable_gaussian_slab()`](https://vandenman.github.io/PDMPSamplersR/reference/dependent_slab_priors.md)
  [`loglinear_gaussian_scale_slab()`](https://vandenman.github.io/PDMPSamplersR/reference/dependent_slab_priors.md)
  : Dependent slab prior constructors

## brms interface

- [`brm_pdmp()`](https://vandenman.github.io/PDMPSamplersR/reference/brm_pdmp.md)
  : Fit a brms model using PDMP samplers

## Variable selection

- [`model_averaged_mean()`](https://vandenman.github.io/PDMPSamplersR/reference/model_averaged_mean.md)
  : Model-averaged posterior means for stickable coefficients
- [`inclusion_prob()`](https://vandenman.github.io/PDMPSamplersR/reference/inclusion_prob.md)
  : Coefficient-wise posterior inclusion probabilities
- [`check_inclusion_stability()`](https://vandenman.github.io/PDMPSamplersR/reference/check_inclusion_stability.md)
  : Per-chain stability check for inclusion probability estimates
- [`inclusion_bf()`](https://vandenman.github.io/PDMPSamplersR/reference/inclusion_bf.md)
  : Inclusion Bayes factors for stickable coefficients
- [`plot_inclusion()`](https://vandenman.github.io/PDMPSamplersR/reference/plot_inclusion.md)
  : Plot posterior inclusion probabilities
- [`median_probability_model()`](https://vandenman.github.io/PDMPSamplersR/reference/median_probability_model.md)
  : Median probability model

## Setup

- [`pdmpsamplers_setup()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmpsamplers_setup.md)
  : Set up the Julia environment for PDMPSamplers.jl
- [`pdmpsamplers_update()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmpsamplers_update.md)
  : Update the Julia project for PDMPSamplers.jl
- [`pdmpsamplers_status()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmpsamplers_status.md)
  : Show the status of the PDMPSamplers.jl installation
- [`use_local_pdmpsamplers()`](https://vandenman.github.io/PDMPSamplersR/reference/use_local_pdmpsamplers.md)
  : Use a local version of PDMPSamplers.jl
- [`check_for_julia_setup()`](https://vandenman.github.io/PDMPSamplersR/reference/check_for_julia_setup.md)
  : Check if the Julia Project is Setup Properly
- [`cleanup_worker_state()`](https://vandenman.github.io/PDMPSamplersR/reference/cleanup_worker_state.md)
  : Release sampling state in a reused worker
