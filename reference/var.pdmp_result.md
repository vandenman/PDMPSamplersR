# Continuous-time variance of PDMP trace

When `chain = NULL`, the pooled variance across chains is computed as
the mean of within-chain variances plus the between-chain variance of
the chain means.

## Usage

``` r
# S3 method for class 'pdmp_result'
var(x, transforms = NULL, chain = NULL, ...)
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
