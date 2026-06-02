# Support Boundary Control

Build a control list for support-boundary handling in PDMP samplers.

## Usage

``` r
support_boundary_control(
  mode = c("error", "line_search", "line_search_truncated_refresh"),
  max_bisection_steps = 60L,
  time_rtol = 1e-08,
  time_atol = 1e-10,
  clip_fraction = 1 - 1e-10,
  max_refresh_attempts = 20L,
  refresh_probe_time = 1e-04,
  min_safe_time = 1e-12
)
```

## Arguments

- mode:

  Character string, how to handle support-boundary violations where the
  target or gradient becomes undefined during forward trajectory
  probing. One of `"error"` (default, fail fast), `"line_search"`
  (localize the first invalid time via bisection, then error), or
  `"line_search_truncated_refresh"` (heuristic BPS-family recovery that
  first searches for ordinary events before the localized boundary,
  handles such an event if one occurs first, and otherwise refreshes
  velocity at a valid interior point).

- max_bisection_steps:

  Integer, maximum bisection iterations.

- time_rtol:

  Numeric, relative tolerance for bisection.

- time_atol:

  Numeric, absolute tolerance for bisection.

- clip_fraction:

  Numeric in `(0, 1]`, fraction of the last-valid time used as a safe
  interior point after localization. The truncated-refresh mode may
  apply an additional conservative cap before refreshing.

- max_refresh_attempts:

  Integer, maximum number of refreshed velocities to try if
  `mode = "line_search_truncated_refresh"` reaches the support boundary
  before an ordinary event.

- refresh_probe_time:

  Numeric, short forward probe time used to reject immediately invalid
  refreshed velocities. A value of zero disables the probe; otherwise
  the Julia fallback may use a slightly longer scale-aware probe.

- min_safe_time:

  Numeric, minimum time gap used when clipping away from the localized
  boundary.

## Value

A list suitable for the `support_boundary` argument.

## Details

Both `"line_search"` and `"line_search_truncated_refresh"` probe along
the linear ray \\x_0 + t v\\ and are therefore only valid for
BPS/ZigZag-family flows with linear dynamics. For non-linear flows
(e.g., Boomerang) these modes fall back to `"error"` behavior.
