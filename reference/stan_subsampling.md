# Compile and configure subsampling custom-Stan models

\`pdmp_subsample_hpp_path()\` returns the installed thread-local
subset-hook header. \`compile_pdmp_stan_model()\` compiles a Stan source
with the required \`–allow-undefined\` and \`USER_HEADER\` settings.

Specifies exact subsampling for an ordinal Markov random field (OMRF)
likelihood with one likelihood contribution per person. The residual of
a person's gradient is taken relative to an analytic, undamped Hessian
control variate at the anchor, and every grid cell gets a certified roof
on the summed residual rates. Threshold and interaction arguments are
Stan unconstrained block names and are resolved at sampler
initialization. The prior gradient is analytic, so \`prior_standata\`
passed to \[stan_subsampling()\] must follow one of the package's OMRF
prior contracts.

## Usage

``` r
pdmp_subsample_hpp_path()

compile_pdmp_stan_model(path, compile_options = stan_compile_options())

stan_residual_envelope(weights, growth_rates = NULL)

stan_subsampling(
  n_observations,
  subsample_size,
  prior_standata,
  residual_envelope,
  anchor = NULL
)

omrf_residual_envelope(
  X,
  seen,
  thresholds,
  interactions,
  exact_selected_refinement = TRUE,
  person_signed_roof = FALSE,
  signed_only_person_roof = FALSE,
  joint_signed_group_size = 0L,
  shared_current_categorical = FALSE
)

stan_parameter_mapping(path_to_stanmodel, standata, parameters = NULL)
```

## Arguments

- path:

  Path to a Stan source file.

- compile_options:

  Settings from
  [`stan_compile_options`](https://vandenman.github.io/PDMPSamplersR/reference/stan_compile_options.md).

- weights:

  Nonnegative component-by-observation curvature weights.

- growth_rates:

  Optional nonnegative growth coefficient per component.

- n_observations:

  Total number of subsampling likelihood contributions.

- subsample_size:

  Number drawn without replacement at each proposal.

- prior_standata:

  Prior-only data list or JSON path for the same model.

- residual_envelope:

  A certified object from \[stan_residual_envelope()\] or a recognized
  package provider such as \[omrf_residual_envelope()\].

- anchor:

  Optional unconstrained anchor.

- X:

  Integer person-by-node matrix encoded from zero.

- seen:

  Number of observed categories per node.

- thresholds:

  Threshold parameter block name.

- interactions:

  Interaction parameter block name; edge order is the upper triangle
  \`(1,2), (1,3), ..., (P-1,P)\`.

- exact_selected_refinement:

  If \`TRUE\` (default), the selected person's residual is evaluated
  once while thinning against the subset bound and reused for the
  candidate rate. \`FALSE\` thins against the drawn roof only.

- person_signed_roof:

  Linear BPS certificate policy. Signed node residual endpoints are
  summed within each person before taking an absolute bound; the
  factor-summed certificate remains a valid fallback. The estimator and
  retained target are unchanged. Other flows keep the ordinary roof.

- signed_only_person_roof:

  Linear BPS policy which builds the certified signed whole-person
  endpoint roof without also building the factor-summed fallback table.
  It preserves the estimator and uses pointwise factor refinement for
  the selected person. Requires \`person_signed_roof = TRUE\`.

- joint_signed_group_size:

  Linear BPS joint signed marked roof. Use \`4\` for fixed seeded groups
  of four, \`16\`, \`32\` or \`64\` for weighted ragged groups, or \`0\`
  to keep the ordinary scheduler. Requires \`signed_only_person_roof =
  TRUE\`.

- shared_current_categorical:

  Evaluator that reuses current categorical moments for exact
  active-neighbour patterns. Requires joint groups.

- path_to_stanmodel:

  Stan source or compiled-library path.

- standata:

  Named data list or JSON path.

- parameters:

  Optional unconstrained names, block prefixes, or indices.

## Value

The header path, or the compiled shared-library path.

A declarative residual-envelope specification.

A validated observation-subsampling specification.

A certified Stan residual-envelope specification.

A data frame with resolved 1-based indices and names.
