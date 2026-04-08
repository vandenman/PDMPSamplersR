# Posterior expected values for subsampled custom families

Post-processing functions used by \[brms::posterior_epred()\] for
subsampled custom families. Each function returns the posterior expected
value.

## Usage

``` r
posterior_epred_subsampled_bernoulli_logit(prep)

posterior_epred_subsampled_poisson_log(prep)

posterior_epred_subsampled_gaussian(prep)

posterior_epred_subsampled_negbinomial_log(prep)

posterior_epred_subsampled_student(prep)

posterior_epred_subsampled_gamma(prep)

posterior_epred_subsampled_beta(prep)
```

## Arguments

- prep:

  A preparation object as used by brms post-processing functions.
