# Online summaries of the retained trajectory

Configures exact continuous-time summaries that are accumulated while a
streaming chain runs. Instead of storing every event, the sampler keeps,
for \`batches\` equal physical-time batches of the retained phase, the
integrals of \\x\\, \\x^2\\, \\x^4\\ and of the free (non-zero)
indicator for every coordinate, the integrals of the number of free
coordinates and its square, and the number of sticky transitions per
coordinate. Optionally it also records the state on a regular time grid.
Memory and disk use therefore do not grow with the number of events.

## Usage

``` r
online_summary_options(
  batches = 400L,
  grid_spacing = NA_real_,
  write_events = FALSE,
  warmup_batches = 0L
)
```

## Arguments

- batches:

  Number of equal physical-time batches of the retained phase.

- grid_spacing:

  Physical time between state snapshots, or \`NA\` for none.

- write_events:

  If \`FALSE\` (default), no retained-phase events are written; the
  online summaries are the only retained output.

- warmup_batches:

  Number of equal physical-time batches of the warmup phase with the
  same integrals (no snapshots), or 0 (default) for none. They are read
  with \`read_online_summaries(path, phase = "warmup")\`.

## Value

A list of class \`pdmp_online_summaries\`.

## Details

Pass the result as \`online_summaries\` to \[pdmp_sample()\] or
\[pdmp_sample_from_stanmodel()\] together with \`trace_storage =
"streaming"\`, and read the summaries afterwards with
\[read_online_summaries()\].
