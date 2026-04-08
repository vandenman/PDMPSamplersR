# Show the Stan code used by brm_pdmp

Returns the Stan code that \[brm_pdmp()\] would compile, after any
subsampling modifications. Without \`subsample_size\`, this is identical
to \`brms::stancode()\`. With \`subsample_size\`, two variants are
returned:

- standard:

  The standard model (used for the full-data and prior-only models).

- ext_cpp:

  The model with likelihood rewritten to use external C++ subsetting
  functions (used for the subsampled-gradient model).

## Usage

``` r
brm_stancode(
  formula,
  data,
  family = gaussian(),
  prior = NULL,
  subsample_size = NULL,
  stanvars = NULL,
  sample_prior = "no",
  ...
)
```

## Arguments

- formula:

  A brms model formula.

- data:

  A data frame containing the variables in the model.

- family:

  A family object (e.g., \`gaussian()\`, \`bernoulli()\`).

- prior:

  A \`brmsprior\` object or NULL for default priors.

- subsample_size:

  Integer number of observations per subsample, or NULL (default) for
  full-data gradients. When non-NULL, a BridgeStan control-variate
  subsampled gradient is used. Must be less than \`nrow(data)\`.
  Fixed-effects and random-effects models are supported for subsampling.

- stanvars:

  Optional \`stanvar\` object for custom Stan code.

- sample_prior:

  Currently only \`"no"\` is supported.

- ...:

  Additional arguments passed to \[brms::brm()\] for model setup.

## Value

A character string when \`subsample_size\` is NULL, or a named list with
elements \`standard\` and \`ext_cpp\`.
