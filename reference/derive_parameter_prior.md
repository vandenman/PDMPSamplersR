# Derive parameter_prior from brms prior specification

Evaluates the slab density at zero for each unconstrained coordinate.
Only \`normal(0, s)\` and \`student_t(df, 0, s)\` priors on \`b\` class
coefficients are supported. All other cases error out.

## Usage

``` r
derive_parameter_prior(
  prior,
  unc_names,
  can_stick,
  supported_coef_names = NULL
)
```

## Arguments

- prior:

  A \`brmsprior\` data frame (from \`brms::prior_summary()\` or the
  prior slot of a brmsfit).

- unc_names:

  Character vector from \`BridgeStan::param_unc_names()\`.

- can_stick:

  Logical vector of length \`d\` (from \`map_can_stick()\`).

- supported_coef_names:

  Optional character vector of supported coefficient names from brms
  metadata.

## Value

Numeric vector of length \`d\` with slab densities at zero for stickable
coordinates and 1.0 elsewhere.
