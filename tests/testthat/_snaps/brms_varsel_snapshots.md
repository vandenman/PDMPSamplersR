# snapshot: can_stick for gaussian y ~ x1 + x2 + x3

    Code
      data.frame(parameter = unc_names, can_stick = result)
    Output
          parameter can_stick
      1 b.Intercept     FALSE
      2        b.x1      TRUE
      3        b.x2      TRUE
      4        b.x3      TRUE
      5       sigma     FALSE

# snapshot: can_stick for logistic y ~ x1 + x2

    Code
      data.frame(parameter = unc_names, can_stick = result)
    Output
          parameter can_stick
      1 b.Intercept     FALSE
      2        b.x1      TRUE
      3        b.x2      TRUE

# snapshot: can_stick for gaussian with scale parameter

    Code
      data.frame(parameter = unc_names, can_stick = result)
    Output
          parameter can_stick
      1 b.Intercept     FALSE
      2        b.x1      TRUE
      3        b.x2      TRUE
      4        b.x3      TRUE
      5        b.x4      TRUE
      6        b.x5      TRUE
      7       sigma     FALSE

# snapshot: stickable_coef_names for 10-predictor model

    Code
      result
    Output
       [1] "b.x1"  "b.x2"  "b.x3"  "b.x4"  "b.x5"  "b.x6"  "b.x7"  "b.x8"  "b.x9" 
      [10] "b.x10"

# snapshot: can_stick with user override

    Code
      data.frame(parameter = unc_names, can_stick = result)
    Output
          parameter can_stick
      1 b.Intercept     FALSE
      2        b.x1      TRUE
      3        b.x2     FALSE
      4        b.x3      TRUE
      5       sigma     FALSE

# snapshot: metadata-first mapping for simple numeric formula

    Code
      data.frame(parameter = unc_names, can_stick = result)
    Output
          parameter can_stick
      1 b.Intercept     FALSE
      2        b.x1      TRUE
      3        b.x2      TRUE
      4        b.f1     FALSE
      5       sigma     FALSE

# snapshot: metadata-first mapping rejects factor formula

    Code
      PDMPSamplersR:::supported_b_coef_names(formula = y ~ f, data = df, fe_names = c(
        "f_b"))
    Condition
      Error in `PDMPSamplersR:::supported_b_coef_names()`:
      ! Formula-based sticky auto-mapping supports numeric predictors only.
      x Non-numeric terms: "f".
      i Factors and grouped structures are not supported.

