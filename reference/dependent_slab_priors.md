# Dependent slab prior constructors

These constructors describe slab priors for the dependent aggregate
sticky Julia path. They are intentionally distinct from the legacy
`kappa` / `parameter_prior` sticky interface: `slab_prior` uses the new
boundary-density convention in which outgoing boundary velocities are
drawn by the aggregate sticky kernel.

## Usage

``` r
dense_gaussian_slab(mean, cov, coef = NULL)

exchangeable_gaussian_slab(mean = 0, u, v, coef = NULL)

loglinear_gaussian_scale_slab(
  log_base_sd,
  logscale,
  logscale_design,
  coef = NULL
)
```

## Arguments

- mean:

  Common mean for an exchangeable Gaussian slab.

- cov:

  Positive-definite covariance matrix for a fixed Gaussian slab.

- coef:

  Optional character vector of unconstrained coefficient names, or
  integer vector of 1-based unconstrained parameter indices, controlled
  by this slab. If omitted, the bridge uses `can_stick`.

- u:

  Diagonal covariance component in \\u I + v 11'\\.

- v:

  Rank-one covariance component in \\u I + v 11'\\.

- log_base_sd:

  Numeric log base standard deviations. A scalar is recycled over the
  selected coefficients.

- logscale:

  Integer or character vector of unconstrained log-scale coordinate
  indices/names.

- logscale_design:

  Numeric coefficient-by-logscale design matrix. Sparse matrices are
  accepted and retained as a compact specification.

## Value

A slab prior specification for the dependent aggregate sticky bridge.
