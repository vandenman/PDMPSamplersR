# Continuous-time covariance of PDMP trace

Continuous-time covariance of PDMP trace

## Usage

``` r
# S3 method for class 'pdmp_result'
cov(x, chain = 1L, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- chain:

  Integer, which chain to use (default: 1).

- ...:

  Ignored.

## Value

Numeric matrix of size `x$d` by `x$d`.
