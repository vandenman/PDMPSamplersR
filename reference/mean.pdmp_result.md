# Continuous-time mean of PDMP trace

Computes the time-weighted mean from the continuous occupation measure.
When `transforms` is supplied, returns the mean of constrained
(transformed) parameters.

## Usage

``` r
# S3 method for class 'pdmp_result'
mean(x, transforms = NULL, chain = NULL, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- transforms:

  Optional list of transform specifications (see
  [`identity_transform`](https://vandenman.github.io/PDMPSamplersR/reference/identity_transform.md)).

- chain:

  Integer, which chain to use, or `NULL` to pool across all chains
  (default: `NULL`).

- ...:

  Ignored.

## Value

Numeric vector of length `x$d`.
