# Plot posterior inclusion probabilities

Produces a horizontal bar chart of coefficient-wise inclusion
probabilities from a sticky \`brmsfit\`, with a vertical reference line
at \`threshold\`.

## Usage

``` r
plot_inclusion(x, ...)

# S3 method for class 'brmsfit'
plot_inclusion(x, threshold = 0.5, ...)
```

## Arguments

- x:

  A \`brmsfit\` fitted with \`sticky = TRUE\`.

- ...:

  Passed to \[inclusion_prob()\].

- threshold:

  Numeric scalar in \\\[0, 1\]\\, the reference line position (default:
  0.5).

## Value

A \`ggplot\` object.
