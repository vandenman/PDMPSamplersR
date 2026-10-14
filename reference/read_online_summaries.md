# Read online summaries written by a streaming chain

Read online summaries written by a streaming chain

## Usage

``` r
read_online_summaries(path, phase = c("main", "warmup"))
```

## Arguments

- path:

  Either the \`online_summaries.bin\` file or a streaming trace
  directory, which is searched for exactly one such file of \`phase\`.

- phase:

  \`"main"\` (default) for the retained phase or \`"warmup"\`; only used
  when \`path\` is a directory.

## Value

A list with the batch integrals (\`x\`, \`x2\`, \`x4\`, \`free\`: one
column per batch), \`batch_time\`, \`size\`, \`size2\` (integrals of the
number of free coordinates and its square), \`transitions\` (per
coordinate and batch), the snapshot grid (\`grid_times\`, \`grid_x\`,
\`grid_free\`), and \`start_time\`, \`width\` and \`end_time\` of the
retained phase.
