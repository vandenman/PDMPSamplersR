# Show the Stan data used by brm_pdmp

Returns the Stan data that \[brm_pdmp()\] would pass to BridgeStan,
after any subsampling modifications. Without \`subsample_size\`, this is
identical to \`brms::standata()\`. With \`subsample_size\`, returns a
named list with three elements:

- full:

  Full dataset.

- prior:

  Prior-only dummy data (\`N=1\`, \`prior_only=1\`).

- subsample:

  An initial subsample of \`subsample_size\` rows.

## Usage

``` r
brm_standata(
  formula,
  data,
  family = gaussian(),
  prior = NULL,
  subsample_size = NULL,
  indices = NULL,
  stanvars = NULL,
  sample_prior = "no",
  ...
)
```

## Arguments

- formula:

  A brms model formula.

- data:

  A data frame containing the variables in the model.

- family:

  A family object (e.g., \`gaussian()\`, \`bernoulli()\`).

- prior:

  A \`brmsprior\` object or NULL for default priors.

- subsample_size:

  Integer number of observations per subsample, or NULL (default) for
  full-data gradients. When non-NULL, a BridgeStan control-variate
  subsampled gradient is used. Must be less than \`nrow(data)\`.
  Fixed-effects and random-effects models are supported for subsampling.

- indices:

  Optional integer vector of observation indices for the initial
  subsample. If NULL (default), a random sample of size
  \`subsample_size\` is drawn.

- stanvars:

  Optional \`stanvar\` object for custom Stan code.

- sample_prior:

  Currently only \`"no"\` is supported.

- ...:

  Additional arguments passed to \[brms::brm()\] for model setup.

## Value

A list (Stan data) when \`subsample_size\` is NULL, or a named list with
elements \`full\`, \`prior\`, and \`subsample\`.
