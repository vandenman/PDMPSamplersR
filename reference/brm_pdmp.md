# Fit a brms model using PDMP samplers

Uses PDMPSamplers.jl as a sampling backend for brms models. Returns a
standard \`brmsfit\` object with all post-processing (\`summary\`,
\`plot\`, \`conditional_effects\`, \`loo\`, etc.) working.

## Usage

``` r
brm_pdmp(
  formula,
  data,
  family = gaussian(),
  prior = NULL,
  flow = c("ZigZag", "BouncyParticle", "Boomerang", "AdaptiveBoomerang",
    "PreconditionedZigZag", "PreconditionedBPS", "DensePreconditionedZigZag",
    "DensePreconditionedBPS"),
  algorithm = c("GridThinningStrategy", "ThinningStrategy", "RootsPoissonStrategy"),
  adaptive_scheme = c("diagonal", "fullrank", "lowrank"),
  T = 50000,
  t0 = 0,
  t_warmup = 0,
  flow_mean = NULL,
  flow_cov = NULL,
  c0 = 0.01,
  grid_n = 30,
  grid_t_max = 2,
  show_progress = TRUE,
  discretize_dt = NULL,
  n_chains = 1L,
  threaded = FALSE,
  seed = NULL,
  compute_lp = FALSE,
  subsample_size = NULL,
  n_anchor_updates = 0L,
  use_hcv = FALSE,
  use_anchor_bank = FALSE,
  bank_capacity = 20L,
  curvature_backend = c("exact", "finite_difference"),
  sticky = FALSE,
  can_stick = NULL,
  model_prior = NULL,
  kappa = NULL,
  slab_prior = NULL,
  stanvars = NULL,
  sample_prior = "no",
  save_model = NULL,
  control = sampler_control(),
  compile_options = stan_compile_options(),
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

- flow:

  Character string specifying the PDMP flow type.

- algorithm:

  Character string specifying the Poisson time strategy.

- adaptive_scheme:

  Character string for AdaptiveBoomerang covariance adaptation:
  \`"diagonal"\`, \`"fullrank"\`, or \`"lowrank"\`.

- T:

  Numeric total simulation time.

- t0:

  Numeric start time.

- t_warmup:

  Numeric warmup duration. Auto-set for adaptive flows. When subsampling
  is active and \`t_warmup\` is 0, it is automatically set to 20% of the
  sampling time.

- flow_mean:

  Numeric vector for the flow reference mean, or NULL. For observation
  subsampling this is also the fixed control-variate anchor and initial
  position; supplying a posterior mode or other representative point can
  materially tighten the residual envelope.

- flow_cov:

  Numeric matrix for the flow covariance, or NULL.

- c0:

  Numeric thinning bound constant. With \`ThinningStrategy\` this
  user-configured bound must dominate the deterministic subsampling
  component; violations are detected and stop the sampler.

- grid_n:

  Integer number of grid points for GridThinningStrategy.

- grid_t_max:

  Numeric max grid interval for GridThinningStrategy.

- show_progress:

  Logical; show sampling progress bar.

- discretize_dt:

  Numeric time step for discretization, or NULL for automatic (yields
  ~1000 samples).

- n_chains:

  Integer number of chains (default: 1). With multiple chains, \`Rhat\`
  and multi-chain diagnostics become available.

- threaded:

  Logical; run chains in parallel (default: FALSE).

- seed:

  NULL (default) or a non-negative integer seed passed through to
  Julia's sampler RNG.

- compute_lp:

  Logical; compute \`lp\_\_\` via \`BridgeStan::log_density()\` for each
  sample (default: FALSE). Adds overhead but enables
  \`bridge_sampler()\` and populates the \`lp\_\_\` diagnostic column.

- subsample_size:

  Integer number of observations per subsample, or NULL (default) for
  full-data gradients. With \`GridThinningStrategy\` or a representable
  \`ThinningStrategy\` envelope, subsampling acceleration is selected
  automatically when affine predictor, first-batch likelihood, and
  trajectory providers are all available for the requested dynamics.
  Ineligible models use exact full gradients.

- n_anchor_updates:

  Integer number of anchor updates during warmup (default: 0). For
  analytic subsampling sampling, each update prepares a coherent anchor
  state containing the full and deterministic anchor gradients, cached
  likelihood predictors, and residual envelope.

- use_hcv:

  Logical; enable the certified analytic damped Hessian control variate
  for subsampling affine Bernoulli-logit or binomial-logit sampling
  (default: FALSE). Its likelihood Hessian and Taylor-remainder envelope
  are prepared with each anchor, with Cauchy damping \`d / (d + \|\|x -
  anchor\|\|^2)\` in \`d\` unconstrained dimensions; finite-difference
  HCV is not used.

- use_anchor_bank:

  Logical; retain prepared warmup anchors and select the nearest one at
  event boundaries (default: FALSE). Selection is chain-local and
  performs no full-data gradient or envelope rebuild. Requires
  \`n_anchor_updates \> 0\`; use only \`n_anchor_updates\` for the
  lower-memory single-anchor update mode.

- bank_capacity:

  Integer capacity of the anchor bank (default: 20). Least-recently-used
  entries are replaced during warmup. Each populated entry stores
  predictor and kernel caches of order \`N \* K\` plus a complete
  residual envelope of order \`N \* R\`; memory therefore grows linearly
  in the populated capacity and can be substantial for large \`N\` or
  many growth-envelope components. Benchmark representative models
  before using a large capacity. Analytic HCV additionally stores one
  dense \`d \* d\` likelihood Hessian per populated anchor.

- curvature_backend:

  Directional curvature evaluation: \`"exact"\` (default) uses
  Hessian-vector products; \`"finite_difference"\` uses shifted
  gradients and has finite-difference truncation error.

- sticky:

  Logical; enable spike-and-slab variable selection for population-level
  coefficients (default: FALSE). Requires \`model_prior\` to be set.
  Sticky transitions redraw the complete velocity from the exact law on
  each active coordinate stratum; no pre-freeze velocity is restored.
  This applies to every supported dynamics, including
  dense-preconditioned Zig-Zag.

- can_stick:

  Optional logical vector indicating which non-intercept
  population-level coefficients are candidates for selection. Length
  must match the number of supported coefficients. If omitted, all
  non-intercept population-level coefficients are candidates.

- model_prior:

  A \[bernoulli()\] or \[betabernoulli()\] object specifying the prior
  on model space.

- kappa:

  Optional numeric vector of slab densities at zero for each stickable
  coordinate (kappa in the sticky PDMP literature). If omitted, derived
  automatically from the brms prior specification (only \`normal(0, s)\`
  and \`student_t(df, 0, s)\` are supported for automatic derivation).

- slab_prior:

  Optional dependent slab prior created by \[dense_gaussian_slab()\],
  \[exchangeable_gaussian_slab()\], or
  \[loglinear_gaussian_scale_slab()\]. Mutually exclusive with
  \`kappa\`. Full-data dependent slabs use the brms prior-only data as
  the base prior target; subsampled brms dependent slabs are not yet
  supported. Adaptive Boomerang is supported with
  \`GridThinningStrategy\` and requires positive warmup (automatically
  set to one fifth of the sampling time when omitted).

- stanvars:

  Optional \`stanvar\` object for custom Stan code.

- sample_prior:

  Currently only \`"no"\` is supported.

- save_model:

  Optional file path to save the generated Stan code.

- control:

  Additional call-local settings from
  [`sampler_control`](https://vandenman.github.io/PDMPSamplersR/reference/sampler_control.md).

- compile_options:

  Settings from
  [`stan_compile_options`](https://vandenman.github.io/PDMPSamplersR/reference/stan_compile_options.md).

- ...:

  Additional arguments passed to \[brms::brm()\] for model setup.

## Value

A \`brmsfit\` object.

## Details

PDMP samplers do not produce NUTS-style diagnostics. The \`lp\_\_\`,
\`accept_stat\_\_\`, and related diagnostic columns are set to zero.
Functions relying on NUTS diagnostics (e.g., \`pairs()\` divergence
plots, \`nuts_params()\`) will not produce meaningful output.

\`loo()\` works because brms computes \`log_lik\` directly from model
parameters, not from \`lp\_\_\`.

With a single chain, \`Rhat\` will report \`NA\`. Use \`n_chains \>= 2\`
for convergence diagnostics.

Subsampling acceleration is enabled when the model has certified affine
predictor, likelihood, and dynamics providers. The first provider batch
covers weighted and subsetted Bernoulli/binomial logit models,
categorical/multinomial logit, Poisson log models (including
\`rate()\`), and Gaussian location-scale models (including known
\`se()\` values and affine sigma predictors). Independent affine
multivariate Bernoulli responses are also supported. Unknown geometry
falls back to the exact full-gradient sampler. Priors, Jacobians, and
unconditional custom target additions remain in an opaque deterministic
BridgeStan provider.

PDMPSamplers.jl applies \`N / m\`, draws a fresh subsampling subset per
proposal, uses dynamics-specific linear or harmonic trajectory bounds,
and reuses the accepted stochastic gradient for the event. Observation
subsampling supports \`GridThinningStrategy\` generally and
\`ThinningStrategy\` when the residual envelope has an affine global
roof (or the trajectory is periodic).
