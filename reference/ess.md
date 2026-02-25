# Effective sample size

Estimate ESS per coordinate from a PDMP trace using the batch-means
method, without discretization.

## Usage

``` r
ess(x, ...)

# S3 method for class 'pdmp_result'
ess(x, chain = 1L, n_batches = 0L, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- ...:

  Passed to methods.

- chain:

  Integer, which chain to use (default: 1).

- n_batches:

  Integer, number of batches. 0 uses the default (default: 0).

## Value

Numeric vector of length `x$d`.
