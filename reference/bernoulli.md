# Create a Bernoulli model prior

Specifies independent Bernoulli inclusion probabilities for use as a
model prior in spike-and-slab sampling.

## Usage

``` r
bernoulli(prob = 0.5)
```

## Arguments

- prob:

  Numeric value or vector of inclusion probabilities, each between 0 and
  1 (default: 0.5). A single value is recycled across all dimensions.

## Value

An object of class `"bernoulli"`.

## See also

[`betabernoulli`](https://vandenman.github.io/PDMPSamplersR/reference/betabernoulli.md)
for a Beta-Bernoulli alternative.
