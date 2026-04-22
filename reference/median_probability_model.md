# Median probability model

Returns the names of population-level coefficients with posterior
inclusion probability above \`threshold\` (default 0.5). This is the
\*median probability model\*: the model that includes each variable for
which the posterior probability of inclusion exceeds one half.

## Usage

``` r
median_probability_model(x, ...)

# S3 method for class 'brmsfit'
median_probability_model(x, threshold = 0.5, ...)
```

## Arguments

- x:

  A \`brmsfit\` fitted with \`sticky = TRUE\`.

- ...:

  Passed to \[inclusion_prob()\].

- threshold:

  Numeric scalar in \\\[0, 1\]\\ (default: 0.5).

## Value

A character vector of coefficient names.

## Collinearity caveat

The median probability model can behave poorly under moderate or strong
collinearity. Coefficient-wise inclusion and term-wise predictive
reduction are not the same task; see the package documentation for
guidance.
