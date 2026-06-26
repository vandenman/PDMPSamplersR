# Empirical CDF

Compute the fraction of trajectory time with the transform of coordinate
`coordinate` at or below `q`.

## Usage

``` r
cdf(x, ...)

# S3 method for class 'pdmp_result'
cdf(x, q, coordinate, transforms = NULL, chain = NULL, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- ...:

  Passed to methods.

- q:

  Numeric threshold.

- coordinate:

  Integer coordinate index.

- transforms:

  Optional list of transform specifications.

- chain:

  Integer, which chain to use, or `NULL` to average across all chains
  (default: `NULL`).

## Value

Numeric scalar in \[0, 1\].
