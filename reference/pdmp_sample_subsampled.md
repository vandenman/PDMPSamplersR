# PDMP sampling with subsampled gradients

Performs PDMP sampling using a user-supplied subsampled gradient
callback. This is a low-level interface intended for correctness
verification and prototyping. For large-scale problems with cheap
likelihoods, consider the Julia-native built-in targets (when available)
to avoid R-to-Julia call overhead.

## Usage

``` r
pdmp_sample_subsampled(
  grad_sub,
  n_obs,
  d,
  subsample_size,
  flow = c("ZigZag", "BouncyParticle", "Boomerang", "AdaptiveBoomerang",
    "PreconditionedZigZag", "PreconditionedBPS"),
  algorithm = c("GridThinningStrategy", "ThinningStrategy", "RootsPoissonStrategy"),
  T = 50000,
  t0 = 0,
  t_warmup = 0,
  flow_mean = NULL,
  flow_cov = NULL,
  c0 = 0.01,
  x0 = NULL,
  hvp_sub = NULL,
  grad_full = NULL,
  use_full_gradient_for_reflections = FALSE,
  grid_n = 30,
  grid_t_max = 2,
  show_progress = TRUE,
  n_chains = 1L,
  threaded = FALSE,
  adaptive_scheme = c("diagonal", "fullrank")
)
```

## Arguments

- grad_sub:

  Function that computes the NEGATIVE subsampled gradient. Should take
  two arguments: a numeric vector \`x\` (position) and an integer vector
  \`indices\` (1-based observation indices), and return a numeric vector
  of length \`d\`.

- n_obs:

  Integer, total number of observations in the dataset.

- d:

  Integer, dimension of the parameter space.

- subsample_size:

  Integer, number of observations per subsample.

- flow:

  Character string specifying the flow type. One of \`"ZigZag"\`,
  \`"BouncyParticle"\`, \`"Boomerang"\`, \`"AdaptiveBoomerang"\`,
  \`"PreconditionedZigZag"\`, or \`"PreconditionedBPS"\`.

- algorithm:

  Character string specifying the algorithm. One of
  \`"ThinningStrategy"\`, \`"GridThinningStrategy"\`, or
  \`"RootsPoissonStrategy"\`.

- T:

  Numeric, total sampling time (default: 50000).

- t0:

  Numeric, initial time (default: 0.0).

- t_warmup:

  Numeric, warmup time (default: 0.0).

- flow_mean:

  Numeric vector of length d, mean for the flow (default: zero vector).

- flow_cov:

  Numeric matrix of size d x d, covariance for the flow (default:
  identity).

- c0:

  Numeric, bound parameter (default: 1e-2).

- x0:

  Numeric vector of length d, initial position (default: random normal).

- hvp_sub:

  Function for the subsampled Hessian-vector product (optional). Takes
  \`(x, v, indices)\` and returns a numeric vector of length \`d\`. If
  \`NULL\`, a finite-difference approximation is used.

- grad_full:

  Function for the full (non-subsampled) NEGATIVE gradient (optional).
  Takes a numeric vector \`x\` and returns a numeric vector of length
  \`d\`. Required when \`use_full_gradient_for_reflections = TRUE\`.

- use_full_gradient_for_reflections:

  Logical, whether to use the full gradient for reflection events
  (default: \`FALSE\`).

- grid_n:

  Integer, number of grid points for GridThinningStrategy (default: 30).

- grid_t_max:

  Numeric, maximum time for grid (default: 2.0).

- show_progress:

  Logical, whether to show progress bar (default: \`TRUE\`).

- n_chains:

  Integer, number of chains to run (default: 1).

- threaded:

  Logical, whether to run chains in parallel (default: \`FALSE\`).

- adaptive_scheme:

  Character string, adaptation scheme for AdaptiveBoomerang (default:
  \`"diagonal"\`).

## Value

A `pdmp_result` object.
