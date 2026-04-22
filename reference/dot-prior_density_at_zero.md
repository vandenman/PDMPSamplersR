# Evaluate the slab density at zero for a single b coefficient

Looks up the prior for \`coef\` first as a coefficient-specific prior,
then falls back to the class-level \`b\` prior. Supports \`normal(0,
s)\` and \`student_t(df, 0, s)\` only.

## Usage

``` r
.prior_density_at_zero(prior, coef)
```

## Arguments

- prior:

  A \`brmsprior\` data frame.

- coef:

  Character, the brms coefficient name (without \`b.\` prefix).

## Value

Numeric scalar, the slab density evaluated at zero.
