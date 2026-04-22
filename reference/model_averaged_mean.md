# Model-averaged posterior means for stickable coefficients

Returns \\E\[\beta_k \mid y\]\\ for each stickable population-level
coefficient by averaging the posterior draws produced by the sticky PDMP
trajectory. Because the trajectory spends time at zero (the spike) in
proportion to the posterior spike mass, this equals the spike-and-slab
model-averaged mean provided the sampler has converged and the
zero-dwell time is correctly attributed. Use \[inclusion_prob()\] to
inspect the corresponding inclusion probabilities.

## Usage

``` r
model_averaged_mean(x, ...)

# S3 method for class 'brmsfit'
model_averaged_mean(x, ...)
```

## Arguments

- x:

  A \`brmsfit\` fitted with \`sticky = TRUE\`.

- ...:

  Ignored.

## Value

A named numeric vector of posterior means.
