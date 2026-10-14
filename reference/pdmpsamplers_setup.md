# Set up the Julia environment for PDMPSamplers.jl

Initializes the Julia runtime via JuliaCall and sets up the Julia
project that PDMPSamplersR depends on. This installs all required Julia
packages and precompiles them. Call this once before using any sampling
functions.

## Usage

``` r
pdmpsamplers_setup(verbose = TRUE, project_dir = NULL)
```

## Arguments

- verbose:

  Logical, if `TRUE` (default), progress messages are printed.

- project_dir:

  Optional Julia project directory. Sets the R option
  `pdmpsamplersr_config_julia_project_dir` for this R session; the same
  option can be set before loading the package. Otherwise the user-cache
  project is used. Configure each independent worker session separately.

## Value

Invisible `TRUE` if the setup succeeds.

## See also

[`pdmpsamplers_update`](https://vandenman.github.io/PDMPSamplersR/reference/pdmpsamplers_update.md)
to update Julia dependencies.
