# Posterior predictions for subsampled custom families

Post-processing functions used by \[brms::posterior_predict()\] for
subsampled custom families. Each function draws posterior predictions
for a single observation.

## Usage

``` r
posterior_predict_subsampled_bernoulli_logit(i, prep, ...)

posterior_predict_subsampled_poisson_log(i, prep, ...)

posterior_predict_subsampled_gaussian(i, prep, ...)

posterior_predict_subsampled_negbinomial_log(i, prep, ...)

posterior_predict_subsampled_student(i, prep, ...)

posterior_predict_subsampled_gamma(i, prep, ...)

posterior_predict_subsampled_beta(i, prep, ...)
```

## Arguments

- i:

  Integer index of the observation.

- prep:

  A preparation object as used by brms post-processing functions.

- ...:

  Additional arguments (unused).
