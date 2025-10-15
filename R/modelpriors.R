#'@export
bernoulli <- function(prob = .5) {
    if (any(!is.numeric(prob)) || any(prob < 0) || any(prob > 1))
      cli::cli_abort("Argument 'prob' must be a numeric value between 0 and 1.")

    structure(list(prob = prob), class = "bernoulli")
}

is.bernoulli <- function(x) {
  inherits(x, "bernoulli")
}

#'@export
betabernoulli <- function(a = 1, b = 1) {
    if (a < 0 || b < 0)
      cli::cli_abort("Arguments 'a' and 'b' must be positive numeric values.")

    structure(list(a = a, b = b), class = "beta-bernoulli")
}

is.betabernoulli <- function(x) {
  inherits(x, "beta-bernoulli")
}
