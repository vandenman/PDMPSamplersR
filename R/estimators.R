# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

.ensure_chains <- function(x) {
  if (!is.null(x$chains)) return(x$chains)
  if (is.null(x$skeleton))
    cli::cli_abort(
      c("No live Julia reference and no skeleton found.",
        "i" = "Call {.fn materialize} before {.fn saveRDS} to enable reload.")
    )
  check_for_julia_setup()
  if (isTRUE(x$skeleton[[1L]]$sparse)) {
    .pdmpsamplers_julia_call("r_from_sparse_skeleton",
      lapply(x$skeleton, `[[`, "initial_time"),
      lapply(x$skeleton, `[[`, "initial_position"),
      lapply(x$skeleton, `[[`, "initial_velocity"),
      lapply(x$skeleton, `[[`, "event_indices"),
      lapply(x$skeleton, `[[`, "event_times"),
      lapply(x$skeleton, `[[`, "event_positions"),
      lapply(x$skeleton, `[[`, "event_velocities"),
      lapply(x$skeleton, `[[`, "is_boomerang"),
      lapply(x$skeleton, `[[`, "mu")
    )
  } else {
    .pdmpsamplers_julia_call("r_from_skeleton",
      lapply(x$skeleton, `[[`, "times"),
      lapply(x$skeleton, `[[`, "positions"),
      lapply(x$skeleton, `[[`, "velocities"),
      lapply(x$skeleton, `[[`, "is_boomerang"),
      lapply(x$skeleton, `[[`, "mu"),
      lapply(x$skeleton, `[[`, "is_mutable")
    )
  }
}

.check_pdmp <- function(x) {
  if (!inherits(x, "pdmp_result"))
    cli::cli_abort("{.arg x} must be a {.cls pdmp_result} object.")
}

.transform_specs <- function(transforms) {
  if (is.null(transforms)) return(NULL)
  lapply(transforms, function(tr) {
    if (!is.list(tr) || is.null(tr$type))
      cli::cli_abort("Each transform must be a list with a {.field type} field.")
    tr
  })
}

# ──────────────────────────────────────────────────────────────────────────────
# mean (base generic)
# ──────────────────────────────────────────────────────────────────────────────

#' Continuous-time mean of PDMP trace
#'
#' Computes the time-weighted mean from the continuous occupation measure.
#' When \code{transforms} is supplied, returns the mean of constrained
#' (transformed) parameters.
#'
#' @param x A \code{pdmp_result} object.
#' @param transforms Optional list of transform specifications (see
#'   \code{\link{identity_transform}}).
#' @param chain Integer, which chain to use, or \code{NULL} to pool across
#'   all chains (default: \code{NULL}).
#' @param ... Ignored.
#'
#' @returns Numeric vector of length \code{x$d}.
#' @export
mean.pdmp_result <- function(x, transforms = NULL, chain = NULL, ...) {
  .check_pdmp(x)
  chains <- .ensure_chains(x)
  specs  <- .transform_specs(transforms)

  if (is.null(chain)) {
    n <- length(chains)
    result <- Reduce(`+`, lapply(seq_len(n), function(ch) {
      if (is.null(specs))
        .pdmpsamplers_julia_call("r_mean", chains, chain = ch)
      else
        .pdmpsamplers_julia_call("r_mean", chains, specs, chain = ch)
    })) / n
    return(result)
  }

  chain <- as.integer(chain)
  if (is.null(specs))
    .pdmpsamplers_julia_call("r_mean", chains, chain = chain)
  else
    .pdmpsamplers_julia_call("r_mean", chains, specs, chain = chain)
}

# ──────────────────────────────────────────────────────────────────────────────
# var / sd / cov / cor — custom S3 generics with base fallback
# ──────────────────────────────────────────────────────────────────────────────

#' Variance
#'
#' Generic function for variance computation. Falls back to
#' \code{\link[stats]{var}} for non-PDMP objects.
#'
#' @param x Object.
#' @param ... Passed to methods.
#'
#' @export
var <- function(x, ...) UseMethod("var")

#' @export
var.default <- function(x, ...) stats::var(x, ...)

#' Continuous-time variance of PDMP trace
#'
#' When \code{chain = NULL}, the pooled variance across chains is computed
#' as the mean of within-chain variances plus the between-chain variance
#' of the chain means.
#'
#' @inheritParams mean.pdmp_result
#' @returns Numeric vector of length \code{x$d}.
#' @export
var.pdmp_result <- function(x, transforms = NULL, chain = NULL, ...) {
  .check_pdmp(x)
  chains <- .ensure_chains(x)
  specs  <- .transform_specs(transforms)

  if (is.null(chain)) {
    n <- length(chains)
    chain_vars <- lapply(seq_len(n), function(ch) {
      if (is.null(specs))
        .pdmpsamplers_julia_call("r_var", chains, chain = ch)
      else
        .pdmpsamplers_julia_call("r_var", chains, specs, chain = ch)
    })
    mean_of_vars <- Reduce(`+`, chain_vars) / n
    if (n > 1L) {
      chain_means <- lapply(seq_len(n), function(ch) {
        if (is.null(specs))
          .pdmpsamplers_julia_call("r_mean", chains, chain = ch)
        else
          .pdmpsamplers_julia_call("r_mean", chains, specs, chain = ch)
      })
      mean_of_vars <- mean_of_vars + apply(do.call(rbind, chain_means), 2L, stats::var)
    }
    return(mean_of_vars)
  }

  chain <- as.integer(chain)
  if (is.null(specs))
    .pdmpsamplers_julia_call("r_var", chains, chain = chain)
  else
    .pdmpsamplers_julia_call("r_var", chains, specs, chain = chain)
}

#' Standard deviation
#'
#' @inheritParams var
#' @export
sd <- function(x, ...) UseMethod("sd")

#' @export
sd.default <- function(x, ...) stats::sd(x, ...)

#' @inheritParams mean.pdmp_result
#' @returns Numeric vector of length \code{x$d}.
#' @export
sd.pdmp_result <- function(x, transforms = NULL, chain = NULL, ...) {
  sqrt(var.pdmp_result(x, transforms = transforms, chain = chain, ...))
}

#' Covariance matrix
#'
#' @inheritParams var
#' @export
cov <- function(x, ...) UseMethod("cov")

#' @export
cov.default <- function(x, ...) stats::cov(x, ...)

#' Continuous-time covariance of PDMP trace
#'
#' @param x A \code{pdmp_result} object.
#' @param chain Integer, which chain to use (default: 1).
#' @param ... Ignored.
#'
#' @returns Numeric matrix of size \code{x$d} by \code{x$d}.
#' @export
cov.pdmp_result <- function(x, chain = 1L, ...) {
  .check_pdmp(x)
  .pdmpsamplers_julia_call("r_cov", .ensure_chains(x), chain = as.integer(chain))
}

#' Correlation matrix
#'
#' @inheritParams var
#' @export
cor <- function(x, ...) UseMethod("cor")

#' @export
cor.default <- function(x, ...) stats::cor(x, ...)

#' Continuous-time correlation of PDMP trace
#'
#' @inheritParams cov.pdmp_result
#' @returns Numeric matrix of size \code{x$d} by \code{x$d}.
#' @export
cor.pdmp_result <- function(x, chain = 1L, ...) {
  .check_pdmp(x)
  .pdmpsamplers_julia_call("r_cor", .ensure_chains(x), chain = as.integer(chain))
}

# ──────────────────────────────────────────────────────────────────────────────
# quantile / median (base generics)
# ──────────────────────────────────────────────────────────────────────────────

#' Continuous-time quantile of PDMP trace
#'
#' When \code{chain = NULL}, per-chain quantiles are averaged pointwise
#' across chains.
#'
#' @param x A \code{pdmp_result} object.
#' @param probs Numeric vector of probabilities in (0, 1).
#' @param transforms Optional list of transform specifications.
#' @param chain Integer, which chain to use, or \code{NULL} to pool across
#'   all chains (default: \code{NULL}).
#' @param coordinate Integer coordinate index, or -1 for all (default: -1).
#' @param ... Ignored.
#'
#' @returns Numeric matrix (\code{length(probs)} x d) when \code{coordinate = -1},
#'   or a numeric vector of length \code{length(probs)} for a single coordinate.
#' @importFrom stats quantile
#' @export
quantile.pdmp_result <- function(x, probs, transforms = NULL, chain = NULL, coordinate = -1L, ...) {
  .check_pdmp(x)
  chains <- .ensure_chains(x)
  specs  <- .transform_specs(transforms)
  coordinate <- as.integer(coordinate)

  if (is.null(chain)) {
    n <- length(chains)
    result_list <- lapply(seq_len(n), function(ch) {
      quantile.pdmp_result(x, probs = probs, transforms = transforms,
                           chain = as.integer(ch), coordinate = coordinate)
    })
    pooled <- Reduce(`+`, result_list) / n
    if (is.matrix(pooled)) rownames(pooled) <- NULL
    return(pooled)
  }

  chain <- as.integer(chain)

  # Julia r_quantile only accepts scalar probs; loop for vectors
  if (length(probs) > 1L) {
    result_list <- lapply(probs, function(p) {
      if (is.null(specs))
        .pdmpsamplers_julia_call("r_quantile", chains, p, chain = chain, coordinate = coordinate)
      else
        .pdmpsamplers_julia_call("r_quantile", chains, p, specs, chain = chain, coordinate = coordinate)
    })
    result <- do.call(rbind, result_list)
    rownames(result) <- NULL
    return(result)
  }

  if (is.null(specs))
    .pdmpsamplers_julia_call("r_quantile", chains, probs, chain = chain, coordinate = coordinate)
  else
    .pdmpsamplers_julia_call("r_quantile", chains, probs, specs, chain = chain, coordinate = coordinate)
}

#' Continuous-time median of PDMP trace
#'
#' @inheritParams quantile.pdmp_result
#' @param na.rm Ignored (present for compatibility with the \link[stats]{median}).
#' @returns Numeric vector or scalar.
#' @importFrom stats median
#' @export
median.pdmp_result <- function(x, na.rm = FALSE, transforms = NULL, chain = NULL, coordinate = -1L, ...) {
  .check_pdmp(x)
  chains <- .ensure_chains(x)
  specs  <- .transform_specs(transforms)
  coordinate <- as.integer(coordinate)

  if (is.null(chain)) {
    n <- length(chains)
    result <- Reduce(`+`, lapply(seq_len(n), function(ch) {
      median.pdmp_result(x, transforms = transforms, chain = as.integer(ch), coordinate = coordinate)
    })) / n
    return(result)
  }

  chain <- as.integer(chain)
  if (is.null(specs))
    .pdmpsamplers_julia_call("r_median", chains, chain = chain, coordinate = coordinate)
  else
    .pdmpsamplers_julia_call("r_median", chains, specs, chain = chain, coordinate = coordinate)
}

# ──────────────────────────────────────────────────────────────────────────────
# Non-base functions: ess, cdf, inclusion_probs, discretize
# ──────────────────────────────────────────────────────────────────────────────

#' Effective sample size
#'
#' Estimate ESS per coordinate from a PDMP trace using the batch-means
#' method, without discretization.
#'
#' @param x A \code{pdmp_result} object.
#' @param ... Passed to methods.
#'
#' @returns Numeric vector of length \code{x$d}.
#' @export
ess <- function(x, ...) UseMethod("ess")

#' @rdname ess
#' @param chain Integer, which chain to use, or \code{NULL} to sum across
#'   all chains (default: \code{NULL}).
#' @param n_batches Integer, number of batches. 0 uses the default (default: 0).
#' @export
ess.pdmp_result <- function(x, chain = NULL, n_batches = 0L, ...) {
  .check_pdmp(x)
  chains <- .ensure_chains(x)

  if (is.null(chain)) {
    n <- length(chains)
    result <- Reduce(`+`, lapply(seq_len(n), function(ch) {
      .pdmpsamplers_julia_call("r_ess", chains, chain = ch, n_batches = as.integer(n_batches))
    }))
    return(result)
  }

  .pdmpsamplers_julia_call("r_ess", chains, chain = as.integer(chain), n_batches = as.integer(n_batches))
}

#' Empirical CDF
#'
#' Compute the fraction of trajectory time with the transform of
#' coordinate \code{coordinate} at or below \code{q}.
#'
#' @param x A \code{pdmp_result} object.
#' @param ... Passed to methods.
#'
#' @returns Numeric scalar in [0, 1].
#' @export
cdf <- function(x, ...) UseMethod("cdf")

#' @rdname cdf
#' @param q Numeric threshold.
#' @param coordinate Integer coordinate index.
#' @param transforms Optional list of transform specifications.
#' @param chain Integer, which chain to use, or \code{NULL} to average
#'   across all chains (default: \code{NULL}).
#' @export
cdf.pdmp_result <- function(x, q, coordinate, transforms = NULL, chain = NULL, ...) {
  .check_pdmp(x)
  chains <- .ensure_chains(x)
  specs  <- .transform_specs(transforms)
  coordinate <- as.integer(coordinate)

  if (is.null(chain)) {
    n <- length(chains)
    result <- Reduce(`+`, lapply(seq_len(n), function(ch) {
      cdf.pdmp_result(x, q = q, coordinate = coordinate, transforms = transforms,
                      chain = as.integer(ch))
    })) / n
    return(result)
  }

  chain <- as.integer(chain)
  if (is.null(specs))
    .pdmpsamplers_julia_call("r_cdf", chains, q, chain = chain, coordinate = coordinate)
  else
    .pdmpsamplers_julia_call("r_cdf", chains, q, specs, chain = chain, coordinate = coordinate)
}

#' Marginal inclusion probabilities
#'
#' For spike-and-slab models, compute the fraction of time each
#' coordinate spent away from zero.
#'
#' @param x A \code{pdmp_result} object.
#' @param ... Passed to methods.
#'
#' @returns Numeric vector of length \code{x$d}.
#' @export
inclusion_probs <- function(x, ...) UseMethod("inclusion_probs")

#' @rdname inclusion_probs
#' @param chain Integer, which chain to use, or \code{NULL} to average
#'   across all chains (default: \code{NULL}).
#' @export
inclusion_probs.pdmp_result <- function(x, chain = NULL, ...) {
  .check_pdmp(x)
  chains <- .ensure_chains(x)

  if (is.null(chain)) {
    n <- length(chains)
    result <- Reduce(`+`, lapply(seq_len(n), function(ch) {
      .pdmpsamplers_julia_call("r_inclusion_probs", chains, chain = ch)
    })) / n
    return(result)
  }

  .pdmpsamplers_julia_call("r_inclusion_probs", chains, chain = as.integer(chain))
}

#' Discretize a PDMP trace
#'
#' Convert the continuous-time trace to a matrix of equally-spaced samples.
#' When \code{dt = NULL} (default), the number of discretization points is set
#' to the continuous-time ESS, preserving the full information content of the
#' trace.
#'
#' @param x A \code{pdmp_result} object.
#' @param ... Passed to methods.
#'
#' @returns Numeric matrix (samples x dimensions).
#' @export
discretize <- function(x, ...) UseMethod("discretize")

#' @rdname discretize
#' @param dt Numeric time step, or NULL to use adaptive discretization based on
#'   continuous-time ESS (default: NULL).
#' @param chain Integer, which chain to use (default: 1).
#' @export
discretize.pdmp_result <- function(x, dt = NULL, chain = 1L, ...) {
  .check_pdmp(x)
  chains <- .ensure_chains(x)
  chain <- as.integer(chain)
  if (is.null(dt))
    .pdmpsamplers_julia_call("r_discretize", chains, chain = chain)
  else
    .pdmpsamplers_julia_call("r_discretize", chains, dt = dt, chain = chain)
}
