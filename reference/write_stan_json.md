# Write data to a JSON file readable by Stan

Write data to a JSON file readable by Stan

## Usage

``` r
write_stan_json(data, file, always_decimal = FALSE)
```

## Arguments

- data:

  (list) A named list of R objects.

- file:

  (string) The path to where the data file should be written.

- always_decimal:

  (logical) Force generate non-integers with decimal points to better
  distinguish between integers and floating point values. If \`TRUE\`
  all R objects in \`data\` intended for integers must be of integer
  type.

## Details

Note: this method is copied from the cmdstanr package to avoid a
dependency on a non-CRAN package. It's copied at commit
\[edccf2d2f6449e7d80626a3ee6cc93845e82915b\](https://github.com/stan-dev/cmdstanr/blob/edccf2d2f6449e7d80626a3ee6cc93845e82915b/R/data.R#L59).
As such, this file and the code therein follows the same license
cmdstanr (BSD 3-Clause License), see the LICENSE file of cmdstanr for
details.

\`write_stan_json()\` performs several conversions before writing the
JSON file:

\* \`logical\` -\> \`integer\` (\`TRUE\` -\> \`1\`, \`FALSE\` -\> \`0\`)
\* \`data.frame\` -\> \`matrix\` (via \[data.matrix()\]) \* \`list\` -\>
\`array\` \* \`table\` -\> \`vector\`, \`matrix\`, or \`array\`
(depending on dimensions of table)

The \`list\` to \`array\` conversion is intended to make it easier to
prepare the data for certain Stan declarations involving arrays:

\* \`vector\[J\] v\[K\]\` (or equivalently \`array\[K\] vector\[J\] v \`
as of Stan 2.27) can be constructed in R as a list with \`K\` elements
where each element a vector of length \`J\` \* \`matrix\[I,J\] v\[K\]\`
(or equivalently \`array\[K\] matrix\[I,J\] m \` as of Stan 2.27 ) can
be constructed in R as a list with \`K\` elements where each element an
\`IxJ\` matrix

These can also be passed in from R as arrays instead of lists but the
list option is provided for convenience. Unfortunately for arrays with
more than one dimension, e.g., \`vector\[J\] v\[K,L\]\` (or equivalently
\`array\[K,L\] vector\[J\] v \` as of Stan 2.27) it is not possible to
use an R list and an array must be used instead. For this example the
array in R should have dimensions \`KxLxJ\`.

## Examples

``` r
x <- matrix(rnorm(10), 5, 2)
y <- rpois(nrow(x), lambda = 10)
z <- c(TRUE, FALSE)
data <- list(N = nrow(x), K = ncol(x), x = x, y = y, z = z)

# write data to json file
file <- tempfile(fileext = ".json")
write_stan_json(data, file)

# check the contents of the file
cat(readLines(file), sep = "\n")
#> {
#>   "N": 5,
#>   "K": 2,
#>   "x": [
#>     [0.25531705484526, -1.82181766097663],
#>     [-2.43726361121953, -0.247325302073524],
#>     [-0.00557128674616073, -0.244199606778383],
#>     [0.621552721415214, -0.282705448814465],
#>     [1.14841160602606, -0.553699383688721]
#>   ],
#>   "y": [11, 16, 4, 11, 7],
#>   "z": [1, 0]
#> }


# demonstrating list to array conversion
# suppose x is declared as `vector[3] x[2]` (or equivalently `array[2] vector[3] x`)
# we can use a list of length 2 where each element is a vector of length 3
data <- list(x = list(1:3, 4:6))
file <- tempfile(fileext = ".json")
write_stan_json(data, file)
cat(readLines(file), sep = "\n")
#> {
#>   "x": [
#>     [1, 2, 3],
#>     [4, 5, 6]
#>   ]
#> }
```
