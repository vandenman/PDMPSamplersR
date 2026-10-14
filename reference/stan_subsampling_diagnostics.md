# Diagnose a subsampling custom-Stan gradient and envelope

Evaluates the persistent full, prior-only, and selected Stan modes at
one position without starting a sampler. The selected likelihood
gradient is unscaled; \`selected_scale\` reports the \`N / m\` factor
applied by Julia when forming the subsampling gradient. The residual
event-rate diagnostic uses the same typed envelope and trajectory
geometry as sampling.

## Usage

``` r
stan_subsampling_diagnostics(
  path_to_stanmodel,
  standata,
  subsampling,
  position,
  subset,
  velocity = NULL,
  flow = c("ZigZag", "BouncyParticle", "AdaptiveBoomerang"),
  flow_mean = NULL,
  flow_cov = NULL,
  phase = c("main", "warmup")
)
```

## Arguments

- path_to_stanmodel:

  Stan source or compiled-library path containing the package subset
  hook.

- standata:

  Full-likelihood Stan data list or JSON path.

- subsampling:

  A specification from \[stan_subsampling()\].

- position:

  Finite unconstrained parameter vector.

- subset:

  One-based observation indices, with exactly the configured subsample
  size and no duplicates.

- velocity:

  Finite diagnostic velocity. Defaults to a vector of ones.

- flow:

  Supported subsampling flow used for the rate diagnostic.

- flow_mean:

  Optional flow mean.

- flow_cov:

  Optional flow covariance.

- phase:

  Diagnostic sampling phase. For staged HCV specifications, \`"main"\`
  activates HCV and \`"warmup"\` uses the first-order fallback.

## Value

A list of gradient-closure, residual-rate, envelope, and construction
diagnostics.
