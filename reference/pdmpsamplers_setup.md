# Set up the Julia environment for PDMPSamplers.jl

Initializes the Julia runtime via JuliaCall and sets up the Julia
project that PDMPSamplersR depends on. This installs all required Julia
packages and precompiles them. Call this once before using any sampling
functions.

## Usage

``` r
pdmpsamplers_setup(verbose = TRUE)
```

## Arguments

- verbose:

  Logical, if `TRUE` (default), progress messages are printed.

## Value

Invisible `TRUE` if the setup succeeds.

## See also

[`pdmpsamplers_update`](https://vandenman.github.io/PDMPSamplersR/reference/pdmpsamplers_update.md)
to update Julia dependencies.
