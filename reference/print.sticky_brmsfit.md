# Print method for sticky brmsfit objects

Delegates to the standard brms print method and appends a one-line note
about sticky variable selection results.

## Usage

``` r
# S3 method for class 'sticky_brmsfit'
print(x, digits = 2, short = getOption("brms.short_summary", FALSE), ...)
```

## Arguments

- x:

  A \`sticky_brmsfit\` object returned by \[brm_pdmp()\] with \`sticky =
  TRUE\`.

- digits:

  Number of decimal digits to print (default: 2).

- short:

  Logical; whether to use the short summary format (default: value of
  \`getOption("brms.short_summary", FALSE)\`).

- ...:

  Passed to the brms print method.

## Value

\`x\`, invisibly.
