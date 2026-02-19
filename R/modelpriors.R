#' Create a Bernoulli model prior
#'
#' Specifies independent Bernoulli inclusion probabilities for use as a
#' model prior in spike-and-slab sampling.
#'
#' @param prob Numeric value or vector of inclusion probabilities, each
#'   between 0 and 1 (default: 0.5). A single value is recycled across
#'   all dimensions.
#'
#' @returns An object of class \code{"bernoulli"}.
#'
#' @seealso \code{\link{betabernoulli}} for a Beta-Bernoulli alternative.
#' @export
bernoulli <- function(prob = .5) {
    if (any(!is.numeric(prob)) || any(prob < 0) || any(prob > 1))
      cli::cli_abort("Argument 'prob' must be a numeric value between 0 and 1.")

    structure(list(prob = prob), class = "bernoulli")
}

is.bernoulli <- function(x) {
  inherits(x, "bernoulli")
}

#' Create a Beta-Bernoulli model prior
#'
#' Specifies a Beta-Bernoulli model prior for spike-and-slab sampling.
#' The inclusion probability is given a \eqn{\text{Beta}(a, b)}{Beta(a, b)}
#' prior, which induces a multiplicity correction.
#'
#' @param a Positive numeric, first shape parameter of the Beta distribution (default: 1).
#' @param b Positive numeric, second shape parameter of the Beta distribution (default: 1).
#'
#' @returns An object of class \code{"beta-bernoulli"}.
#'
#' @seealso \code{\link{bernoulli}} for fixed inclusion probabilities.
#' @export
betabernoulli <- function(a = 1, b = 1) {
    if (a < 0 || b < 0)
      cli::cli_abort("Arguments 'a' and 'b' must be positive numeric values.")

    structure(list(a = a, b = b), class = "beta-bernoulli")
}

is.betabernoulli <- function(x) {
  inherits(x, "beta-bernoulli")
}
