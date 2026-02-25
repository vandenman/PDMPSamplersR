# Discretize a PDMP trace

Convert the continuous-time trace to a matrix of equally-spaced samples.

## Usage

``` r
discretize(x, ...)

# S3 method for class 'pdmp_result'
discretize(x, dt = NULL, chain = 1L, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- ...:

  Passed to methods.

- dt:

  Numeric time step, or NULL to use mean inter-event time (default:
  NULL).

- chain:

  Integer, which chain to use (default: 1).

## Value

Numeric matrix (samples x dimensions).
