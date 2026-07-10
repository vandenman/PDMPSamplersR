#' Dependent slab prior constructors
#'
#' These constructors describe slab priors for the dependent aggregate sticky
#' Julia path. They are intentionally distinct from the legacy \code{kappa} /
#' \code{parameter_prior} sticky interface: \code{slab_prior} uses the new
#' boundary-density convention in which outgoing boundary velocities are drawn
#' by the aggregate sticky kernel.
#'
#' @param kappa Positive slab density at zero. A scalar is recycled by the
#'   Julia bridge over the selected sticky coefficients.
#' @param coef Optional character vector of unconstrained coefficient names, or
#'   integer vector of 1-based unconstrained parameter indices, controlled by
#'   this slab. If omitted, the bridge uses \code{can_stick}.
#'
#' @returns A slab prior specification for the dependent aggregate sticky
#'   bridge.
#'
#' @export
independent_slab_density <- function(kappa, coef = NULL) {
  validate_type(kappa, type = "double", positive = TRUE)
  .new_slab_prior("independent_slab_density", list(kappa = kappa), coef = coef)
}

#' @rdname independent_slab_density
#' @param mean Numeric mean vector for a fixed Gaussian slab.
#' @param cov Positive-definite covariance matrix for a fixed Gaussian slab.
#' @export
dense_gaussian_slab <- function(mean, cov, coef = NULL) {
  validate_type(mean, type = "double")
  validate_type(cov, type = "double", dims = c(length(mean), length(mean)))
  if (!isSymmetric(cov)) {
    cli::cli_abort("Argument {.arg cov} must be symmetric.")
  }
  if (inherits(try(chol(cov), silent = TRUE), "try-error")) {
    cli::cli_abort("Argument {.arg cov} must be positive definite.")
  }
  .new_slab_prior("dense_gaussian", list(mean = mean, cov = cov), coef = coef)
}

#' @rdname independent_slab_density
#' @param mean Common mean for an exchangeable Gaussian slab.
#' @param u Diagonal covariance component in \eqn{u I + v 11'}.
#' @param v Rank-one covariance component in \eqn{u I + v 11'}.
#' @param zero_mean Logical. If \code{TRUE}, use the optimized zero-mean
#'   exchangeable Julia provider.
#' @export
exchangeable_gaussian_slab <- function(mean = 0, u, v, zero_mean = TRUE, coef = NULL) {
  validate_type(mean, type = "double", n = 1)
  validate_type(u, type = "double", n = 1, positive = TRUE)
  validate_type(v, type = "double", n = 1)
  validate_type(zero_mean, type = "logical", n = 1)
  if (isTRUE(zero_mean) && !isTRUE(all.equal(mean, 0))) {
    cli::cli_abort("Argument {.arg mean} must be 0 when {.arg zero_mean} is TRUE.")
  }
  .new_slab_prior(
    "exchangeable_gaussian",
    list(mean = mean, u = u, v = v, zero_mean = zero_mean),
    coef = coef
  )
}

#' @rdname independent_slab_density
#' @param mean_cov R function for state-dependent Gaussian slab callbacks.
#'   It is called as \code{mean_cov(x)} and should return a list with numeric
#'   fields \code{mean} and \code{cov}.
#' @param active_prior_neggrad R function for the active slab
#'   negative-gradient contribution. It is called as
#'   \code{active_prior_neggrad(x, active)}, where \code{active} is a logical
#'   vector over the slab coefficients, and should return a full unconstrained
#'   negative-gradient vector. This callback is required when a callback
#'   Gaussian slab is used for dependent-slab sampling.
#' @export
gaussian_scale_mixture_slab <- function(mean_cov, active_prior_neggrad = NULL, coef = NULL) {
  if (!rlang::is_function(mean_cov)) {
    cli::cli_abort("Argument {.arg mean_cov} must be a function.")
  }
  if (!is.null(active_prior_neggrad) && !rlang::is_function(active_prior_neggrad)) {
    cli::cli_abort("Argument {.arg active_prior_neggrad} must be NULL or a function.")
  }
  .new_slab_prior(
    "callback_gaussian",
    list(mean_cov = mean_cov, active_prior_neggrad = active_prior_neggrad),
    coef = coef
  )
}

#' @rdname independent_slab_density
#' @param log_q_zero R function returning a boundary log-density at zero. It is
#'   called as \code{log_q_zero(x, active, j)}, where \code{j} is a 1-based slab
#'   coefficient index.
#' @export
arbitrary_slab_boundary <- function(log_q_zero, active_prior_neggrad, coef = NULL) {
  if (!rlang::is_function(log_q_zero)) {
    cli::cli_abort("Argument {.arg log_q_zero} must be a function.")
  }
  if (!rlang::is_function(active_prior_neggrad)) {
    cli::cli_abort("Argument {.arg active_prior_neggrad} must be a function.")
  }
  .new_slab_prior(
    "arbitrary_boundary",
    list(log_q_zero = log_q_zero, active_prior_neggrad = active_prior_neggrad),
    coef = coef
  )
}

is.slab_prior <- function(x) {
  inherits(x, "dependent_slab_prior")
}

.new_slab_prior <- function(type, fields, coef = NULL) {
  if (!is.null(coef) && !(is.character(coef) || rlang::is_integerish(coef))) {
    cli::cli_abort("Argument {.arg coef} must be NULL, a character vector, or an integer vector.")
  }
  if (!is.null(coef) && !is.character(coef)) {
    coef <- as.integer(coef)
  }
  if (!is.null(coef) && anyDuplicated(coef)) {
    cli::cli_abort("Argument {.arg coef} cannot contain duplicates.")
  }
  structure(
    c(list(type = type, coef = coef), fields),
    class = c(paste0(type, "_slab_prior"), "dependent_slab_prior")
  )
}
