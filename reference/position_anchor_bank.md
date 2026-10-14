# Frozen position anchor bank from a set of states

Builds the anchor bank used by position-selected OMRF subsampling from
states of a previous (typically full-gradient) chain, for example the
online grid snapshots of \[read_online_summaries()\]. It applies the
same rule the sampler uses for its own warmup: coordinate scales are the
standard deviations of the states (floored at a tenth of their median),
the first anchor is the last state, and further anchors are added by
greedy farthest-point selection in the scaled metric.

## Usage

``` r
position_anchor_bank(
  standata,
  subsampling,
  parameter_names,
  states,
  capacity = 8L,
  metric = c("scaled", "hessian")
)
```

## Arguments

- standata:

  The Stan data of the subsampled fit.

- subsampling:

  The \[stan_subsampling()\] object of the subsampled fit.

- parameter_names:

  Unconstrained parameter names, in sampler order.

- states:

  Numeric matrix with one state per column.

- capacity:

  Maximum number of anchors.

- metric:

  Anchor selection metric, \`"scaled"\` or \`"hessian"\`.

## Value

A version-1 position anchor bank specification.

## Details

During sampling each state uses the anchor nearest in a shared metric,
so the anchor is a function of position and the Voronoi cells are
bounded by hyperplanes. \`metric = "scaled"\` uses the coordinate scales
above; \`metric = "hessian"\` uses the mean likelihood Hessian of the
anchors plus the identity, which approximates the posterior precision
and so predicts the control-variate error better (it requires global HCV
anchors).

Pass the result as \`subsampling_position_anchor_bank_spec\` to
\[pdmp_sample_from_stanmodel()\], with
\`subsampling_anchor_selection_policy = "position"\`, \`t_warmup = 0\`
and a frozen flow.
