# Create a Beta-Bernoulli model prior

Specifies a Beta-Bernoulli model prior for spike-and-slab sampling. The
inclusion probability is given a \\\text{Beta}(a, b)\\ prior, which
induces a multiplicity correction.

## Usage

``` r
betabernoulli(a = 1, b = 1)
```

## Arguments

- a:

  Positive numeric, first shape parameter of the Beta distribution
  (default: 1).

- b:

  Positive numeric, second shape parameter of the Beta distribution
  (default: 1).

## Value

An object of class `"beta-bernoulli"`.

## See also

[`bernoulli`](https://vandenman.github.io/PDMPSamplersR/reference/bernoulli.md)
for fixed inclusion probabilities.
