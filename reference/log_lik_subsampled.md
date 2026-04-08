# Log-likelihood for subsampled custom families

Post-processing functions used by \[brms::log_lik()\] for subsampled
custom families. Each function computes the log-likelihood for a single
observation.

## Usage

``` r
log_lik_subsampled_bernoulli_logit(i, prep)

log_lik_subsampled_poisson_log(i, prep)

log_lik_subsampled_gaussian(i, prep)

log_lik_subsampled_negbinomial_log(i, prep)

log_lik_subsampled_student(i, prep)

log_lik_subsampled_gamma(i, prep)

log_lik_subsampled_beta(i, prep)
```

## Arguments

- i:

  Integer index of the observation.

- prep:

  A preparation object as used by brms post-processing functions.
