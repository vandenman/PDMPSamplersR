# Additional sampler controls

Pass this object as `control` to
[`pdmp_sample()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmp_sample.md),
[`pdmp_sample_from_stanmodel()`](https://vandenman.github.io/PDMPSamplersR/reference/pdmp_sample_from_stanmodel.md)
or
[`brm_pdmp()`](https://vandenman.github.io/PDMPSamplersR/reference/brm_pdmp.md).
Settings apply to that call only, including when a worker is reused for
another model.

## Usage

``` r
sampler_control(
  grid_n_min = NULL,
  grid_bound_violation = NULL,
  adaptive_boomerang_refresh_rate = 0.1,
  warmup_stop = NULL,
  profile = FALSE,
  research_diagnostics = FALSE,
  hcv_cache = TRUE,
  check_fast_deterministic = FALSE
)
```

## Arguments

- grid_n_min:

  Minimum grid size, or `NULL` for `min(grid_n, 5)`.

- grid_bound_violation:

  Grid bound violation policy: `NULL` uses the sampler default;
  otherwise `"throw"`, `"shrink"`, or `"count"`.

- adaptive_boomerang_refresh_rate:

  Initial refresh rate for adaptive Boomerang dynamics. The default is
  0.1.

- warmup_stop:

  Optional named list configuring adaptive warmup stopping: `min_time`,
  `max_time`, `stable_time`, `min_events`, and `check_every`. An empty
  list enables defaults (one quarter, all, and one quarter of
  `t_warmup`, 100 events, and a check every 25 events); `NULL` disables
  it. Supported by custom-Stan sampling.

- profile:

  Collect bridge-phase timings for custom-Stan sampling.

- research_diagnostics:

  Collect additional OMRF allocation and HVP timing diagnostics. These
  measurements add overhead.

- hcv_cache:

  Reuse OMRF Hessian products along a flight.

- check_fast_deterministic:

  Compare the fast OMRF deterministic gradient with its generic
  evaluation. Intended for diagnostics; adds overhead.

## Value

A sampler control list.
