# Marginal inclusion probabilities

For spike-and-slab models, compute the fraction of time each coordinate
spent away from zero.

## Usage

``` r
inclusion_probs(x, ...)

# S3 method for class 'pdmp_result'
inclusion_probs(x, chain = NULL, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- ...:

  Passed to methods.

- chain:

  Integer, which chain to use, or `NULL` to average across all chains
  (default: `NULL`).

## Value

Numeric vector of length `x$d`.
