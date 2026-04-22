# Validate sticky arguments for brm_pdmp

Validate sticky arguments for brm_pdmp

## Usage

``` r
validate_brms_sticky(
  sticky,
  can_stick,
  model_prior,
  parameter_prior,
  d,
  unc_names,
  supported_coef_names = NULL,
  prior,
  subsampled
)
```

## Arguments

- sticky:

  Logical.

- can_stick:

  Logical vector or NULL.

- model_prior:

  A bernoulli/betabernoulli object or NULL.

- parameter_prior:

  Numeric vector or NULL.

- d:

  Integer, total number of unconstrained parameters.

- unc_names:

  Character vector, unconstrained parameter names.

- supported_coef_names:

  Character vector of supported coefficient names derived from brms
  metadata.

- prior:

  brms prior data frame.

- subsampled:

  Logical, whether subsampling is active.

## Value

A list with validated \`can_stick\`, \`model_prior\`, and
\`parameter_prior\` (all length \`d\`).
