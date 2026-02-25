# Check if the Julia Project is Setup Properly

Check if the Julia Project is Setup Properly

## Usage

``` r
check_for_julia_setup(error_if_not_exists = FALSE, setup_if_not_exists = TRUE)
```

## Arguments

- error_if_not_exists:

  logical, if TRUE, will throw an error if the Julia project is not
  found. If FALSE, will attempt to set up the project if
  setup_if_not_exists is TRUE.

- setup_if_not_exists:

  logical, if FALSE and \`error_if_not_exists\` is FALSE, will not throw
  an error but printif the Julia project is not found. If TRUE, will
  attempt to set up the project if \`error_if_not_exists\` is FALSE.

## Value

TRUE if the Julia project is successfully loaded or set up, `FALSE`
otherwise (if \`error_if_not_exists\` is also `FALSE`). If
\`error_if_not_exists\` is `TRUE`, will throw an error if the project
cannot be loaded or set up.

## Details

This function checks for the existence of the Julia project for
PDMPSamplers.jl, and if it exists, loads it and the required packages.
If it does not exist, it can either throw an error or attempt to set up
the project, depending on the parameters. It is mostly useful for
internal use to ensure the Julia environment is ready before calling any
functions that depend on it, but can also be called directly by users if
they want to check or set up the Julia project.
