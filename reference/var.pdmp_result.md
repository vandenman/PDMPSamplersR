# Continuous-time variance of PDMP trace

Continuous-time variance of PDMP trace

## Usage

``` r
# S3 method for class 'pdmp_result'
var(x, transforms = NULL, chain = 1L, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- transforms:

  Optional list of transform specifications (see
  [`identity_transform`](https://vandenman.github.io/PDMPSamplersR/reference/identity_transform.md)).

- chain:

  Integer, which chain to use (default: 1).

- ...:

  Ignored.

## Value

Numeric vector of length `x$d`.
