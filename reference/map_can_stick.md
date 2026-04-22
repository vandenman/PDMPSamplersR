# Build logical can_stick vector from brms metadata

Determines which unconstrained BridgeStan coordinates correspond to
supported population-level coefficients (excluding intercept, scale, and
shape parameters).

## Usage

``` r
map_can_stick(unc_names, supported_coef_names = NULL, user_can_stick = NULL)
```

## Arguments

- unc_names:

  Character vector from \`BridgeStan::param_unc_names()\`.

- supported_coef_names:

  Character vector of supported coefficient names on the unconstrained
  scale (e.g. \`b.x1\`) derived from brms metadata.

- user_can_stick:

  Optional logical vector (same length as supported non-intercept \`b\`
  coefficients) to override the default all-TRUE.

## Value

A logical vector of length \`length(unc_names)\`.
