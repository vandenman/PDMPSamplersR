# Per-chain stability check for inclusion probability estimates

Compares inclusion probability estimates across chains to help assess
convergence of the sticky sampler. Returns a data frame with cross-chain
summary statistics for each stickable coefficient.

## Usage

``` r
check_inclusion_stability(x, ...)

# Default S3 method
check_inclusion_stability(x, ...)

# S3 method for class 'brmsfit'
check_inclusion_stability(x, ...)
```

## Arguments

- x:

  A \`brmsfit\` fitted with \`sticky = TRUE\`.

- ...:

  Ignored.

## Value

A data frame with columns \`coefficient\`, \`mean\`, \`sd\`, \`min\`,
\`max\`, and \`range\` summarising cross-chain variability.
