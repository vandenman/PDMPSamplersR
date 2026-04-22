# Return the names of stickable coefficients

Return the names of stickable coefficients

## Usage

``` r
stickable_coef_names(unc_names, supported_coef_names = NULL)
```

## Arguments

- unc_names:

  Character vector from \`BridgeStan::param_unc_names()\`.

- supported_coef_names:

  Optional character vector of supported coefficient names from brms
  metadata. If omitted, uses legacy raw name matching over
  \`unc_names\`.

## Value

Character vector of stickable coefficient names.
