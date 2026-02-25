# Continuous-time median of PDMP trace

Continuous-time median of PDMP trace

## Usage

``` r
# S3 method for class 'pdmp_result'
median(x, na.rm = FALSE, transforms = NULL, chain = 1L, coordinate = -1L, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- na.rm:

  Ignored (present for compatibility with the
  [median](https://rdrr.io/r/stats/median.html)).

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
