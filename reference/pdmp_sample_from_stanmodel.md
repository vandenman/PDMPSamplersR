# PDMP Sampling from Stan Model

Performs Piecewise Deterministic Markov Process (PDMP) sampling from a
Stan model using the PDMPSamplers.jl Julia package with BridgeStan
integration.

## Usage

``` r
pdmp_sample_from_stanmodel(
  path_to_stanmodel,
  standata,
  subsampling = NULL,
  subsampling_anchor_updates = 0L,
  subsampling_use_anchor_bank = FALSE,
  subsampling_anchor_bank_capacity = 8L,
  subsampling_anchor_selection_policy = c("adaptive", "position"),
  subsampling_position_anchor_bank_spec = NULL,
  subsampling_main_anchor_refresh_distance = Inf,
  subsampling_require_no_main_full_data = FALSE,
  subsampling_diagnostic_window = 0,
  flow = c("ZigZag", "BouncyParticle", "Boomerang", "AdaptiveBoomerang",
    "PreconditionedZigZag", "PreconditionedBPS"),
  algorithm = c("ThinningStrategy", "GridThinningStrategy", "RootsPoissonStrategy"),
  T = 50000,
  t0 = 0,
  t_warmup = 0,
  warmup_adaptation_interval = NULL,
  boomerang_variance_floor_fraction = 0,
  flow_mean = NULL,
  flow_cov = NULL,
  c0 = 0.01,
  flow_specification = c("auto", "initial", "frozen"),
  x0 = NULL,
  theta0 = NULL,
  free0 = NULL,
  stored_velocity0 = NULL,
  sticky = FALSE,
  can_stick = NULL,
  model_prior = NULL,
  parameter_prior = NULL,
  slab_prior = NULL,
  grid_n = 30,
  grid_t_max = 2,
  curvature_backend = NULL,
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
  control = sampler_control(),
  compile_options = stan_compile_options()
)
```

## Arguments

- path_to_stanmodel:

  Character, path to a Stan model file (.stan) or a compiled Stan model
  (.so/.dll/.dylib). If a .stan file is provided, BridgeStan will
  compile it automatically.

- standata:

  Either a character path to the Stan data file (JSON format), or a
  named list that will be written to a temporary JSON file via
  \[write_stan_json()\].

- subsampling:

  Optional exact subsampling custom-Stan specification created by
  \[stan_subsampling()\]. The selected Stan branch must add the unscaled
  sum of the installed observations; Julia applies \`N / m\`.

- subsampling_anchor_updates:

  Number of prepared OMRF anchors added from continuous-time warmup
  means.

- subsampling_use_anchor_bank:

  Retain prepared OMRF anchors and select the nearest entry at event
  boundaries. If false, updates replace the active anchor.

- subsampling_anchor_bank_capacity:

  Maximum prepared OMRF anchors kept per chain.

- subsampling_anchor_selection_policy:

  Research anchor-bank policy. \`"adaptive"\` preserves the existing
  event-history/envelope policy; \`"position"\` selects one frozen
  finite bank source by weighted position and supports linear
  Zig–Zag/BPS and harmonic Boomerang flows.

- subsampling_position_anchor_bank_spec:

  Optional version-1 mathematical position-bank specification saved on a
  previous result. Restoration rebuilds prepared provider quantities
  from its ordered positions and frozen distance scales; it requires
  zero warmup and matching dataset, parameter-order, and
  envelope/grouping identities.

- subsampling_main_anchor_refresh_distance:

  Positive Euclidean distance from the nearest prepared anchor that
  triggers a new exact OMRF anchor during main sampling. The default
  \`Inf\` disables main-phase preparation.

- subsampling_require_no_main_full_data:

  If \`TRUE\`, fail the run when retained sampling performs any
  full-gradient evaluation, prepares an anchor, or traverses every OMRF
  person-node conditional in native Julia. Setup and warmup work are
  allowed. This is intended as a hard benchmark invariant for
  frozen-control-variate and direct-clock experiments.

- subsampling_diagnostic_window:

  Width, in sampling time, of the windows over which the OMRF anchor
  manager records retained-phase subsampling diagnostics. The default 0
  disables them.

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

- warmup_adaptation_interval:

  Nonnegative time between dynamics adaptation attempts during warmup.
  \`NULL\` uses ten equally spaced warmup windows (\`t_warmup / 10\`);
  zero attempts adaptation at every event. The same interval is used by
  full-gradient and subsampling runs.

- boomerang_variance_floor_fraction:

  Nonnegative opt-in research safeguard for adaptive diagonal Boomerang
  flows. At every adaptation update, named \`global_scale\` and
  \`node_scale.\*\` coordinates retain at least this fraction of their
  preadaptation variance. Zero (the default) disables the safeguard.

- flow_mean:

  Numeric vector of length d, mean vector for the flow (default: zero
  vector).

- flow_cov:

  Numeric matrix of size d x d, covariance matrix for the flow (default:
  identity matrix).

- c0:

  Numeric, bound parameter (default: 1e-2).

- flow_specification:

  Role of an explicitly supplied \`flow_mean\` or \`flow_cov\`.
  \`"initial"\` permits legitimate warmup adaptation and records every
  chain's retained flow. \`"frozen"\` requires every retained chain to
  match the request and fails closed otherwise. \`"auto"\` uses
  \`"initial"\` when \`t_warmup \> 0\` and \`"frozen"\` otherwise.

- x0:

  Numeric vector of length d, initial position (default: random normal).

- theta0:

  Optional initial velocity forwarded to both full-gradient and
  custom-Stan subsampling initialization.

- free0, stored_velocity0:

  Optional sticky state to continue a previous run: `free0` is the
  logical mask of free coordinates and `stored_velocity0` the stored
  velocities of the frozen ones. Both must be supplied together, with
  length equal to the number of parameters, and require `sticky = TRUE`.

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

- curvature_backend:

  Character string selecting directional curvature: \`"exact"\` or
  \`"finite_difference"\`. Custom-Stan subsampling has no deterministic
  exact HVP unless using the analytic dependent OMRF provider. Finite
  differences use shifted exact gradients but approximate their
  directional derivative and therefore have step-size and truncation
  error. \`NULL\` selects exact curvature for full gradients and finite
  differences for subsampling, except that analytic dependent OMRF
  always uses its exact composed HVP.

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

- compile_options:

  Settings from
  [`stan_compile_options`](https://vandenman.github.io/PDMPSamplersR/reference/stan_compile_options.md).

## Value

A `pdmp_result` object. Use `mean`, `var`, `quantile`, etc. for
continuous-time estimators, or `discretize` to obtain a sample matrix.

## Details

The gradient and Hessian-vector product are automatically derived from
the Stan model via the BridgeStan extension of PDMPSamplers.jl. This
means you do not need to provide a separate gradient or Hessian
function.
