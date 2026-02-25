# PDMPSamplersR

**PDMPSamplersR** is an R interface to the
[PDMPSamplers.jl](https://github.com/vandenman/PDMPSamplers.jl) Julia
package, enabling Bayesian inference using Piecewise Deterministic
Markov Processes (PDMPs) from within R.

## Status

⚠️ **This package is under active development. Expect breaking changes
and limited support.**

## Installation

`PDMPSamplersR` is only available from GitHub (for now). You can install
it using:

``` r
remotes::install_github("vandenman/PDMPSamplersR")
# or if you use renv:
# renv::install("vandenman/PDMPSamplersR")
```

## Usage

After installation, you can use the R interface to call PDMP-based
samplers via Julia. See the package documentation and vignettes for
examples.

## References

- [PDMPSamplers.jl GitHub
  Repository](https://github.com/vandenman/PDMPSamplers.jl)
