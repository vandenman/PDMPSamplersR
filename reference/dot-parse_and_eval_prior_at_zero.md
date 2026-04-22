# Parse a brms prior string and evaluate the density at zero

Parse a brms prior string and evaluate the density at zero

## Usage

``` r
.parse_and_eval_prior_at_zero(prior_str, coef)
```

## Arguments

- prior_str:

  Character, e.g. "normal(0, 2.5)" or "student_t(3, 0, 2.5)".

- coef:

  Character, coefficient name for error messages.

## Value

Numeric scalar.
