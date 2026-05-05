# Materialise a PDMP result for serialisation

Copies the chain skeleton (event times, positions, velocities) from
Julia into R and drops the live Julia reference. The returned object can
be saved with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and all
estimators continue to work after reloading, provided Julia is set up
before the first estimator call.

## Usage

``` r
materialize(x)
```

## Arguments

- x:

  A `pdmp_result`.

## Value

A `pdmp_result` with `$chains = NULL` and `$skeleton` populated.

## Details

This only needs to be called manually when the result was created with
`materialize = FALSE`; all sampling functions materialise the trace
automatically by default.
