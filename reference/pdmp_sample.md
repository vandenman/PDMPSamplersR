# PDMP Sampling

Performs Piecewise Deterministic Markov Process (PDMP) sampling using
the PDMPSamplers.jl Julia package.

## Usage

``` r
pdmp_sample(
  f,
  d,
  flow = c("ZigZag", "BouncyParticle", "Boomerang"),
  algorithm = c("ThinningStrategy", "GridThinningStrategy", "RootsPoissonStrategy"),
  T = 50000,
  t0 = 0,
  t_warmup = 0,
  flow_mean = NULL,
  flow_cov = NULL,
  c0 = 0.01,
  x0 = NULL,
  theta0 = NULL,
  hessian = NULL,
  sticky = FALSE,
  can_stick = NULL,
  model_prior = NULL,
  parameter_prior = NULL,
  grid_n = 30,
  grid_t_max = 2,
  show_progress = TRUE,
  n_chains = 1L,
  threaded = FALSE
)
```

## Arguments

- f:

  Function that computes the NEGATIVE gradient. Should take a numeric
  vector and return a numeric vector of the same length.

- d:

  Integer, dimension of the problem.

- flow:

  Character string specifying the flow type. One of "ZigZag",
  "BouncyParticle", or "Boomerang".

- algorithm:

  Character string specifying the algorithm. One of "ThinningStrategy",
  "GridThinningStrategy", or "RootsPoissonStrategy".

- T:

  Numeric, total sampling time (default: 50000).

- t0:

  Numeric, initial time (default: 0.0).

- t_warmup:

  Numeric, warmup time (default: 0.0). Events during warmup are
  discarded.

- flow_mean:

  Numeric vector of length d, mean vector for the flow (default: zero
  vector).

- flow_cov:

  Numeric matrix of size d x d, covariance matrix for the flow (default:
  identity matrix).

- c0:

  Numeric, bound parameter (default: 1e-2).

- x0:

  Numeric vector of length d, initial position (default: random normal).

- theta0:

  Numeric vector of length d, initial velocity (default: random based on
  the flow).

- hessian:

  Function that returns the negative Hessian matrix of size d x d
  (default: NULL). Only used with GridThinningStrategy to compute the
  Hessian-vector product.

- sticky:

  Logical, whether to use sticky sampling (default: FALSE).

- can_stick:

  Logical vector of length d, which coordinates can stick (default: all
  FALSE).

- model_prior:

  Prior distribution object for model selection. Should be of class
  'bernoulli' or 'beta-bernoulli' (default: NULL).

- parameter_prior:

  Numeric vector of length d, prior parameters for sticky sampling
  (default: NULL).

- grid_n:

  Integer, number of grid points for GridThinningStrategy (default: 30).

- grid_t_max:

  Numeric, maximum time for grid in GridThinningStrategy (default: 2.0).

- show_progress:

  Logical, whether to show progress bar (default: TRUE).

- n_chains:

  Integer, number of chains to run (default: 1).

- threaded:

  Logical, whether to run chains in parallel (default: FALSE).

## Value

A `pdmp_result` object. Use `mean`, `var`, `quantile`, etc. for
continuous-time estimators, or `discretize` to obtain a sample matrix.
