# Inclusion Bayes factors for stickable coefficients

Computes the inclusion Bayes factor for each stickable coefficient:
\$\$BF\_{\text{incl},k} = \frac{P(\gamma_k = 1 \mid y)}{P(\gamma_k = 0
\mid y)} \cdot \frac{P(\gamma_k = 0)}{P(\gamma_k = 1)}\$\$

## Usage

``` r
inclusion_bf(x, ...)

# Default S3 method
inclusion_bf(x, ...)

# S3 method for class 'brmsfit'
inclusion_bf(x, ...)
```

## Arguments

- x:

  A \`brmsfit\` fitted with \`sticky = TRUE\` using a version of
  \[brm_pdmp()\] that stores the model prior in the sticky metadata.

- ...:

  Passed to \[inclusion_prob()\].

## Value

A named numeric vector of inclusion Bayes factors, one per stickable
coefficient.

## Details

This answers "how much did the data update the prior odds in favour of
inclusion?" Values greater than 1 indicate evidence for inclusion;
values less than 1 indicate evidence against. Unlike inclusion
probabilities alone, the Bayes factor separates the likelihood signal
from the prior.
