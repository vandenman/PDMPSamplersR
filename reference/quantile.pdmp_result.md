# Continuous-time quantile of PDMP trace

Continuous-time quantile of PDMP trace

## Usage

``` r
# S3 method for class 'pdmp_result'
quantile(x, probs, transforms = NULL, chain = 1L, coordinate = -1L, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- probs:

  Numeric scalar in (0, 1).

- transforms:

  Optional list of transform specifications.

- chain:

  Integer, which chain to use (default: 1).

- coordinate:

  Integer coordinate index, or -1 for all (default: -1).

- ...:

  Ignored.

## Value

Numeric vector or scalar.
