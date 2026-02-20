#' Identity transform (no constraint)
#'
#' @returns A transform specification list.
#' @export
identity_transform <- function() {
  list(type = "identity")
}

#' Lower-bounded transform
#'
#' Maps unconstrained \eqn{y \in \mathbb{R}} to \eqn{x = L + \exp(y)}.
#'
#' @param lower Numeric lower bound.
#' @returns A transform specification list.
#' @export
lower_transform <- function(lower) {
  validate_type(lower, type = "double", n = 1)
  list(type = "lower", lower = lower)
}

#' Upper-bounded transform
#'
#' Maps unconstrained \eqn{y \in \mathbb{R}} to \eqn{x = U - \exp(y)}.
#'
#' @param upper Numeric upper bound.
#' @returns A transform specification list.
#' @export
upper_transform <- function(upper) {
  validate_type(upper, type = "double", n = 1)
  list(type = "upper", upper = upper)
}

#' Double-bounded transform
#'
#' Maps unconstrained \eqn{y \in \mathbb{R}} to
#' \eqn{x = L + (U - L) \cdot \mathrm{logistic}(y)}.
#'
#' @param lower Numeric lower bound.
#' @param upper Numeric upper bound.
#' @returns A transform specification list.
#' @export
double_transform <- function(lower, upper) {
  validate_type(lower, type = "double", n = 1)
  validate_type(upper, type = "double", n = 1)
  if (lower >= upper)
    cli::cli_abort("{.arg lower} ({lower}) must be less than {.arg upper} ({upper}).")
  list(type = "double", lower = lower, upper = upper)
}
