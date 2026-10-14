# PDMP Sampling

Performs Piecewise Deterministic Markov Process (PDMP) sampling using
the PDMPSamplers.jl Julia package.

## Usage

``` r
pdmp_sample(
  f,
  d,
  flow = c("ZigZag", "BouncyParticle", "Boomerang", "AdaptiveBoomerang",
    "PreconditionedZigZag", "PreconditionedBPS"),
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
  slab_prior = NULL,
  grid_n = 30,
  grid_t_max = 2,
  grid_bound = c("value_quadratic", "constant"),
  grid_curvature_bound = NULL,
  show_progress = TRUE,
  n_chains = 1L,
  threaded = FALSE,
  seed = NULL,
  adaptive_scheme = c("diagonal", "fullrank"),
  refresh_rate = 1,
  materialize = TRUE,
  trace_storage = c("dense", "streaming"),
  streaming_trace_directory = NULL,
  streaming_trace_buffer_events = 4096L,
  online_summaries = NULL,
  support_boundary = support_boundary_control(),
  control = sampler_control()
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
  "BouncyParticle", "Boomerang", "AdaptiveBoomerang",
  "PreconditionedZigZag", or "PreconditionedBPS". The
  \`"AdaptiveBoomerang"\` flow learns its reference (mean and precision)
  during warmup and requires a grid-like algorithm. It is supported with
  dependent slabs when \`GridThinningStrategy\` is selected; provide
  positive \`t_warmup\`, or let the sampler choose its 20 The
  \`"PreconditionedZigZag"\` and \`"PreconditionedBPS"\` flows learn a
  diagonal preconditioner during warmup and also require a grid-like
  algorithm.

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

  Prior distribution object for model selection:
  [`bernoulli()`](https://vandenman.github.io/PDMPSamplersR/reference/bernoulli.md)
  or
  [`betabernoulli()`](https://vandenman.github.io/PDMPSamplersR/reference/betabernoulli.md).

- parameter_prior:

  Numeric vector of length d, prior parameters for sticky sampling
  (default: NULL).

- slab_prior:

  Optional dependent slab prior created by
  [`dense_gaussian_slab()`](https://vandenman.github.io/PDMPSamplersR/reference/dependent_slab_priors.md),
  [`exchangeable_gaussian_slab()`](https://vandenman.github.io/PDMPSamplersR/reference/dependent_slab_priors.md),
  or
  [`loglinear_gaussian_scale_slab()`](https://vandenman.github.io/PDMPSamplersR/reference/dependent_slab_priors.md).
  Mutually exclusive with `parameter_prior`. The full target supplied
  through `f` must contain the complete coherent slab represented by
  `slab_prior`; the package replaces that slab by its active-face
  restriction internally.

- grid_n:

  Integer, number of grid points for GridThinningStrategy (default: 30).

- grid_t_max:

  Numeric, maximum time for grid in GridThinningStrategy (default: 2.0).

- grid_bound:

  Character string selecting the GridThinning proposal envelope.
  `"value_quadratic"` (default) bounds each grid cell by the larger
  endpoint rate plus `grid_curvature_bound * h^2 / 8`; without
  `grid_curvature_bound` the cells use the ordinary GridThinning
  construction. `"constant"` always uses the ordinary construction.

- grid_curvature_bound:

  Optional nonnegative bound on the curvature of the event rate over a
  grid cell, used by `grid_bound = "value_quadratic"`. Too small a value
  invalidates the envelope.

- show_progress:

  Logical, whether to show progress bar (default: TRUE).

- n_chains:

  Integer, number of chains to run (default: 1).

- threaded:

  Logical, whether to run chains in parallel (default: FALSE).

- seed:

  NULL (default) or a non-negative integer seed. The seed controls
  Julia's sampler RNG and any generated default initial position when
  `x0 = NULL`.

- adaptive_scheme:

  Character string, adaptation scheme for AdaptiveBoomerang. One of
  "diagonal" (default, O(d) per update) or "fullrank" (O(d^3) per
  update, better for correlated targets). Ignored for other flow types.
  With sticky sampling, "fullrank" releases frozen coordinates with the
  experimental boundary-proposal law; "diagonal" uses the
  stored-velocity law.

- refresh_rate:

  Refreshment rate of the `"BouncyParticle"` and `"PreconditionedBPS"`
  flows (default 1). Ignored by the other flows.

- materialize:

  Logical. If `TRUE` (default), the chain skeleton is immediately
  extracted from Julia into R so that
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html) work without any
  extra steps. Set to `FALSE` to skip the extraction and keep only the
  live Julia reference. This can save time and memory if you don't need
  to save the result or if you plan to call
  [`materialize()`](https://vandenman.github.io/PDMPSamplersR/reference/materialize.md)
  manually later.

- trace_storage:

  `"dense"` (default) keeps the event trace in memory. `"streaming"`
  writes it in chunks to `streaming_trace_directory`, so memory use does
  not grow with the number of events; it is required for
  `online_summaries`.

- streaming_trace_directory:

  Directory for the streamed trace. Required when
  `trace_storage = "streaming"`, ignored otherwise.

- streaming_trace_buffer_events:

  Number of events buffered in memory before a chunk is written (default
  4096). Used only with `trace_storage = "streaming"`.

- online_summaries:

  `NULL` (default) or options from
  [`online_summary_options()`](https://vandenman.github.io/PDMPSamplersR/reference/online_summary_options.md).
  Requires `trace_storage = "streaming"`; the summaries are read with
  [`read_online_summaries()`](https://vandenman.github.io/PDMPSamplersR/reference/read_online_summaries.md).

- support_boundary:

  A list created by
  [`support_boundary_control()`](https://vandenman.github.io/PDMPSamplersR/reference/support_boundary_control.md)
  that controls support-boundary diagnostics and heuristic event/refresh
  recovery.

- control:

  Additional call-local settings from
  [`sampler_control`](https://vandenman.github.io/PDMPSamplersR/reference/sampler_control.md).

## Value

A `pdmp_result` object. Use `mean`, `var`, `quantile`, etc. for
continuous-time estimators, or `discretize` to obtain a sample matrix.
