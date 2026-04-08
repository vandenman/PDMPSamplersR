# Discretize a PDMP trace

Convert the continuous-time trace to a matrix of equally-spaced samples.
When `dt = NULL` (default), the number of discretization points is set
to the continuous-time ESS, preserving the full information content of
the trace.

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

  Numeric time step, or NULL to use adaptive discretization based on
  continuous-time ESS (default: NULL).

- chain:

  Integer, which chain to use (default: 1).

## Value

Numeric matrix (samples x dimensions).
