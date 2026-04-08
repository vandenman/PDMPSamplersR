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
    "PreconditionedZigZag", "PreconditionedBPS"),
  algorithm = c("GridThinningStrategy", "ThinningStrategy", "RootsPoissonStrategy"),
  adaptive_scheme = c("diagonal", "fullrank"),
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
  compute_lp = FALSE,
  subsample_size = NULL,
  n_anchor_updates = 10L,
  resample_dt = NULL,
  hvp_mode = c("scaled", "none"),
  use_hcv = FALSE,
  use_anchor_bank = FALSE,
  bank_capacity = 20L,
  use_fd_hvp = FALSE,
  post_warmup_simplify = FALSE,
  use_fd_hcv = FALSE,
  stanvars = NULL,
  sample_prior = "no",
  save_model = NULL,
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

  Character string for preconditioner adaptation.

- T:

  Numeric total simulation time.

- t0:

  Numeric start time.

- t_warmup:

  Numeric warmup duration. Auto-set for adaptive flows. When subsampling
  is active and \`t_warmup\` is 0, it is automatically set to 20% of the
  sampling time.

- flow_mean:

  Numeric vector for the flow reference mean, or NULL.

- flow_cov:

  Numeric matrix for the flow covariance, or NULL.

- c0:

  Numeric thinning bound constant.

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

- compute_lp:

  Logical; compute \`lp\_\_\` via \`BridgeStan::log_density()\` for each
  sample (default: FALSE). Adds overhead but enables
  \`bridge_sampler()\` and populates the \`lp\_\_\` diagnostic column.

- subsample_size:

  Integer number of observations per subsample, or NULL (default) for
  full-data gradients. When non-NULL, a BridgeStan control-variate
  subsampled gradient is used. Must be less than \`nrow(data)\`.
  Fixed-effects and random-effects models are supported for subsampling.

- n_anchor_updates:

  Integer number of anchor updates during warmup (default: 10). Only
  used when \`subsample_size\` is non-NULL.

- resample_dt:

  Numeric time step for resampling, or NULL (default) for no resampling.
  Only used when \`subsample_size\` is non-NULL.

- hvp_mode:

  Character string controlling Hessian-vector product scaling in
  subsampled gradients. One of \`"scaled"\` (default) or \`"none"\`.

- use_hcv:

  Logical; enable damped Hessian control variate (HCV) correction for
  subsampled gradients (default: FALSE). Requires \`subsample_size\` to
  be set. Adds a second-order correction that reduces gradient variance
  near the anchor.

- use_anchor_bank:

  Logical; enable anchor bank with LRU eviction for subsampled gradients
  (default: FALSE). Requires \`subsample_size\` to be set. Maintains
  multiple cached anchor points and selects the nearest one at each
  trajectory step.

- bank_capacity:

  Integer capacity of the anchor bank (default: 20). Only used when
  \`use_anchor_bank\` is TRUE.

- use_fd_hvp:

  Logical; use finite-difference directional curvature instead of
  BridgeStan's Hessian-vector product (default: FALSE). Replaces each
  HVP call (2-3x gradient cost) with one extra gradient call, giving
  ~30-50% per-grid-point savings.

- post_warmup_simplify:

  Logical; switch to a fast constant-bound thinning mode after warmup
  when the sampler is well-adapted (default: FALSE). Activates when the
  reflection ratio is below 30% and at least 10 events have been
  observed.

- use_fd_hcv:

  Logical; use finite-difference approximation for the HVP inside the
  HCV correction (default: FALSE). Replaces the exact BridgeStan HVP at
  the anchor with a cheaper gradient-based FD approximation. Also
  disables HVP for grid thinning. Only used when \`use_hcv\` is TRUE.

- stanvars:

  Optional \`stanvar\` object for custom Stan code.

- sample_prior:

  Currently only \`"no"\` is supported.

- save_model:

  Optional file path to save the generated Stan code.

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

When \`subsample_size\` is specified, the function uses BridgeStan
data-swapping to compute control-variate subsampled gradients. A
centering fix is applied to the brms-generated Stan code so that
predictor centering remains consistent across subsamples.
