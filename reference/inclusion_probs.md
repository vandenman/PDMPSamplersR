# Marginal inclusion probabilities

For spike-and-slab models, compute the fraction of time each coordinate
spent away from zero.

## Usage

``` r
inclusion_probs(x, ...)

# S3 method for class 'pdmp_result'
inclusion_probs(x, chain = 1L, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- ...:

  Passed to methods.

- chain:

  Integer, which chain to use (default: 1).

## Value

Numeric vector of length `x$d`.
