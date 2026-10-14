#' Dependent slab prior constructors
#'
#' These constructors describe slab priors for the dependent aggregate sticky
#' Julia path. They are intentionally distinct from the legacy \code{kappa} /
#' \code{parameter_prior} sticky interface: \code{slab_prior} uses the new
#' boundary-density convention in which outgoing boundary velocities are drawn
#' by the aggregate sticky kernel.
#'
#' @param mean Numeric mean vector for a fixed Gaussian slab.
#' @param cov Positive-definite covariance matrix for a fixed Gaussian slab.
#' @param coef Optional character vector of unconstrained coefficient names, or
#'   integer vector of 1-based unconstrained parameter indices, controlled by
#'   this slab. If omitted, the bridge uses \code{can_stick}.
#'
#' @returns A slab prior specification for the dependent aggregate sticky
#'   bridge.
#'
#' @rdname dependent_slab_priors
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

#' @rdname dependent_slab_priors
#' @param mean Common mean for an exchangeable Gaussian slab.
#' @param u Diagonal covariance component in \eqn{u I + v 11'}.
#' @param v Rank-one covariance component in \eqn{u I + v 11'}.
#' @export
exchangeable_gaussian_slab <- function(mean = 0, u, v, coef = NULL) {
  validate_type(mean, type = "double", n = 1)
  validate_type(u, type = "double", n = 1, positive = TRUE)
  validate_type(v, type = "double", n = 1)
  .new_slab_prior(
    "exchangeable_gaussian",
    list(mean = mean, u = u, v = v),
    coef = coef
  )
}

#' @rdname dependent_slab_priors
#' @param log_base_sd Numeric log base standard deviations. A scalar is
#'   recycled over the selected coefficients.
#' @param logscale Integer or character vector of unconstrained log-scale
#'   coordinate indices/names.
#' @param logscale_design Numeric coefficient-by-logscale design matrix. Sparse
#'   matrices are accepted and retained as a compact specification.
#' @export
loglinear_gaussian_scale_slab <- function(log_base_sd, logscale,
                                           logscale_design, coef = NULL) {
  validate_type(log_base_sd, type = "double")
  if (!(is.character(logscale) || rlang::is_integerish(logscale))) {
    cli::cli_abort("Argument {.arg logscale} must be a character or integer vector.")
  }
  if (!is.matrix(logscale_design) && !inherits(logscale_design, "Matrix")) {
    cli::cli_abort("Argument {.arg logscale_design} must be a numeric matrix or sparse Matrix.")
  }
  design_values <- if (inherits(logscale_design, "Matrix")) {
    Matrix::summary(logscale_design)$x
  } else logscale_design
  if (!is.numeric(design_values) || any(!is.finite(design_values))) {
    cli::cli_abort("Argument {.arg logscale_design} must contain finite numeric values.")
  }
  if (!is.character(logscale) && ncol(logscale_design) != length(logscale)) {
    cli::cli_abort("Columns of {.arg logscale_design} must match {.arg logscale}.")
  }
  if (!is.character(logscale)) logscale <- as.integer(logscale)
  dims <- as.integer(dim(logscale_design))
  design <- if (inherits(logscale_design, "Matrix")) {
    # Adding a declared general sparse zero canonicalizes CSC, triplet,
    # row-sparse, and symmetric Matrix subclasses to their complete matrix
    # semantics before extracting the exported triplet summary.
    general <- Matrix::drop0(logscale_design + Matrix::sparseMatrix(
      i = integer(), j = integer(), dims = dims
    ))
    entries <- Matrix::summary(general)
    list(
      storage = "sparse_rows",
      dims = dims,
      i = as.integer(entries$i),
      j = as.integer(entries$j),
      x = as.numeric(entries$x)
    )
  } else {
    list(storage = "dense", dims = dims, x = as.numeric(logscale_design))
  }
  .new_slab_prior(
    "loglinear_gaussian_scale",
    list(log_base_scales = log_base_sd, logscale = logscale,
         logscale_design = design),
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
