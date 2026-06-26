# Continuous-time quantile of PDMP trace

When `chain = NULL`, per-chain quantiles are averaged pointwise across
chains.

## Usage

``` r
# S3 method for class 'pdmp_result'
quantile(x, probs, transforms = NULL, chain = NULL, coordinate = -1L, ...)
```

## Arguments

- x:

  A `pdmp_result` object.

- probs:

  Numeric vector of probabilities in (0, 1).

- transforms:

  Optional list of transform specifications.

- chain:

  Integer, which chain to use, or `NULL` to pool across all chains
  (default: `NULL`).

- coordinate:

  Integer coordinate index, or -1 for all (default: -1).

- ...:

  Ignored.

## Value

Numeric matrix (`length(probs)` x d) when `coordinate = -1`, or a
numeric vector of length `length(probs)` for a single coordinate.
