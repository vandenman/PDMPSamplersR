#' Write data to a JSON file readable by Stan
#'
#' @export
#' @param data (list) A named list of \R objects.
#' @param file (string) The path to where the data file should be written.
#' @param always_decimal (logical) Force generate non-integers with decimal
#' points to better distinguish between integers and floating point values.
#' If `TRUE` all \R objects in `data` intended for integers must be of integer
#' type.
#'
#' @details
#'
#' Note: this method is copied from the cmdstanr package to avoid a dependency
#' on a non-CRAN package. It's copied at commit [edccf2d2f6449e7d80626a3ee6cc93845e82915b](https://github.com/stan-dev/cmdstanr/blob/edccf2d2f6449e7d80626a3ee6cc93845e82915b/R/data.R#L59).
#' As such, this file and the code therein follows the same license cmdstanr (BSD 3-Clause License),
#' see the LICENSE file of cmdstanr for details.
#'
#' `write_stan_json()` performs several conversions before writing the JSON
#' file:
#'
#' * `logical` -> `integer` (`TRUE` -> `1`, `FALSE` -> `0`)
#' * `data.frame` -> `matrix` (via [data.matrix()])
#' * `list` -> `array`
#' * `table` -> `vector`, `matrix`, or `array` (depending on dimensions of table)
#'
#' The `list` to `array` conversion is intended to make it easier to prepare
#' the data for certain Stan declarations involving arrays:
#'
#' * `vector[J] v[K]` (or equivalently `array[K] vector[J] v ` as of Stan 2.27)
#' can be constructed in \R as a list with `K` elements where each element a
#' vector of length `J`
#' * `matrix[I,J] v[K]` (or equivalently `array[K] matrix[I,J] m ` as of Stan
#' 2.27 ) can be constructed in \R as a list with `K` elements where each element
#' an `IxJ` matrix
#'
#' These can also be passed in from \R as arrays instead of lists but the list
#' option is provided for convenience. Unfortunately for arrays with more than
#' one dimension, e.g., `vector[J] v[K,L]` (or equivalently
#' `array[K,L] vector[J] v ` as of Stan 2.27) it is not possible to use an \R
#' list and an array must be used instead. For this example the array in \R
#' should have dimensions `KxLxJ`.
#'
#' @examples
#' x <- matrix(rnorm(10), 5, 2)
#' y <- rpois(nrow(x), lambda = 10)
#' z <- c(TRUE, FALSE)
#' data <- list(N = nrow(x), K = ncol(x), x = x, y = y, z = z)
#'
#' # write data to json file
#' file <- tempfile(fileext = ".json")
#' write_stan_json(data, file)
#'
#' # check the contents of the file
#' cat(readLines(file), sep = "\n")
#'
#'
#' # demonstrating list to array conversion
#' # suppose x is declared as `vector[3] x[2]` (or equivalently `array[2] vector[3] x`)
#' # we can use a list of length 2 where each element is a vector of length 3
#' data <- list(x = list(1:3, 4:6))
#' file <- tempfile(fileext = ".json")
#' write_stan_json(data, file)
#' cat(readLines(file), sep = "\n")
#'
write_stan_json <- function(data, file, always_decimal = FALSE) {
  if (!is.list(data)) {
    stop("'data' must be a list.", call. = FALSE)
  }
  if (!is.character(file) || !nzchar(file)) {
    stop("The supplied filename is invalid!", call. = FALSE)
  }

  data_names <- names(data)
  if (length(data) > 0 &&
      (length(data_names) == 0 ||
       length(data_names) != sum(nzchar(data_names)))) {
    stop("All elements in 'data' list must have names.", call. = FALSE)

  }
  if (anyDuplicated(data_names) != 0) {
    stop("Duplicate names not allowed in 'data'.", call. = FALSE)
  }

  for (var_name in data_names) {
    var <- data[[var_name]]
    if (!(is.numeric(var) || is.factor(var) || is.logical(var) ||
          is.data.frame(var) || is.list(var))) {
      stop("Variable '", var_name, "' is of invalid type.", call. = FALSE)
    }
    if (anyNA(var)) {
      stop("Variable '", var_name, "' has NA values.", call. = FALSE)
    }

    if (is.table(var)) {
      var <- unclass(var)
    } else if (is.logical(var)) {
      mode(var) <- "integer"
    } else if (is.data.frame(var)) {
      var <- data.matrix(var)
    } else if (is.list(var)) {
      var <- list_to_array(var, var_name)
    }
    data[[var_name]] <- var
  }

  # unboxing variables (N = 10 is stored as N : 10, not N: [10])
  jsonlite::write_json(
    data,
    path = file,
    auto_unbox = TRUE,
    factor = "integer",
    always_decimal = always_decimal,
    digits = NA,
    pretty = TRUE
  )
}


list_to_array <- function(x, name = NULL) {
  list_length <- length(x)
  if (list_length == 0) {
    return(NULL)
  }
  all_dims <- lapply(x, function(z) dim(z) %||% length(z)) # dim is null if vector
  all_equal_dim <- all(sapply(all_dims, function(d) {
    isTRUE(all.equal(d, all_dims[[1]]))
  }))
  if (!all_equal_dim) {
    stop("All matrices/vectors in list '", name, "' must be the same size!", call. = FALSE)
  }
  all_numeric <- all(sapply(x, function(a) is.numeric(a)))
  if (!all_numeric) {
    stop("All elements in list '", name, "' must be numeric!", call. = FALSE)
  }
  element_num_of_dim <- length(all_dims[[1]])
  x <- unlist(x)
  dim(x) <- c(all_dims[[1]], list_length)
  aperm(x, c(element_num_of_dim + 1L, seq_len(element_num_of_dim)))
}
