# Coefficient-wise posterior inclusion probabilities

Extracts per-coefficient inclusion probabilities from a sticky
\`brmsfit\` fitted with \[brm_pdmp()\] and \`sticky = TRUE\`. Each
probability estimates \\P(\gamma_k = 1 \mid y)\\, the posterior
probability that coefficient \\\beta_k\\ is non-zero under the
user-specified spike-and-slab prior.

## Usage

``` r
inclusion_prob(x, ...)

# Default S3 method
inclusion_prob(x, ...)

# S3 method for class 'brmsfit'
inclusion_prob(x, chain = NULL, ...)
```

## Arguments

- x:

  A \`brmsfit\` fitted with \`sticky = TRUE\`.

- ...:

  Ignored.

- chain:

  Integer, which chain to extract, or \`NULL\` to average across all
  chains (default: \`NULL\`).

## Value

A named numeric vector in \\\[0, 1\]\\, one entry per stickable
population-level coefficient.
