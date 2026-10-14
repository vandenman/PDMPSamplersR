# Stan compilation options

Pass as `compile_options` to
[`pdmp_sample_from_stanmodel()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmp_sample_from_stanmodel.md),
[`brm_pdmp()`](https://vandenman.github.io/PDMPSamplersR/reference/brm_pdmp.md)
or
[`compile_pdmp_stan_model()`](https://vandenman.github.io/PDMPSamplersR/reference/stan_subsampling.md).
Existing compiled libraries can be reused; use a separate cache
directory when changing toolchains.

## Usage

``` r
stan_compile_options(
  cache = FALSE,
  cache_dir = NULL,
  native = FALSE,
  optimize = FALSE
)
```

## Arguments

- cache:

  Use a content-addressed compilation cache.

- cache_dir:

  Cache directory, or `NULL` for the Julia depot's
  `pdmpsamplersr/bridgestan` directory.

- native:

  Enable native CPU compiler instructions.

- optimize:

  Enable Stan's `--O1` optimization.

## Value

A compilation options list.
