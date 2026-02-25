# Package index

## Sampling

- [`pdmp_sample()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmp_sample.md)
  : PDMP Sampling
- [`pdmp_sample_from_stanmodel()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmp_sample_from_stanmodel.md)
  : PDMP Sampling from Stan Model
- [`write_stan_json()`](https://vandenman.github.io/PDMPSamplersR/reference/write_stan_json.md)
  : Write data to a JSON file readable by Stan

## Results

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

## Setup

- [`pdmpsamplers_setup()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmpsamplers_setup.md)
  : Set up the Julia environment for PDMPSamplers.jl
- [`pdmpsamplers_update()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmpsamplers_update.md)
  : Update the Julia project for PDMPSamplers.jl
- [`use_local_pdmpsamplers()`](https://vandenman.github.io/PDMPSamplersR/reference/use_local_pdmpsamplers.md)
  : Use a local version of PDMPSamplers.jl
- [`check_for_julia_setup()`](https://vandenman.github.io/PDMPSamplersR/reference/check_for_julia_setup.md)
  : Check if the Julia Project is Setup Properly
