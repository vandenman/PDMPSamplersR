# Continuous-time correlation of PDMP trace

Continuous-time correlation of PDMP trace

## Usage

``` r
# S3 method for class 'pdmp_result'
cor(x, chain = 1L, ...)
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
