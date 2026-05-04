# Constructor (internal)
new_pdmp_result <- function(chains, stats, d, n_chains, skeleton = NULL) {
  structure(
    list(
      chains   = chains,
      stats    = stats,
      d        = d,
      n_chains = n_chains,
      skeleton = skeleton
    ),
    class = "pdmp_result"
  )
}

#' Materialise a PDMP result for serialisation
#'
#' Copies the chain skeleton (event times, positions, velocities) from Julia
#' into R and drops the live Julia reference. The returned object can be saved
#' with \code{saveRDS()} and all estimators continue to work after reloading,
#' provided Julia is set up before the first estimator call.
#'
#' @param x A \code{pdmp_result}.
#' @param ... Ignored.
#' @returns A \code{pdmp_result} with \code{$chains = NULL} and
#'   \code{$skeleton} populated.
#' @export
materialize <- function(x, ...) UseMethod("materialize")

#' @rdname materialize
#' @export
materialize.pdmp_result <- function(x, ...) {
  .check_pdmp(x)
  chains <- if (!is.null(x$chains)) x$chains else JuliaCall::julia_call("r_from_skeleton",
    lapply(x$skeleton, `[[`, "times"),
    lapply(x$skeleton, `[[`, "positions"),
    lapply(x$skeleton, `[[`, "velocities"),
    lapply(x$skeleton, `[[`, "is_boomerang"),
    lapply(x$skeleton, `[[`, "mu"),
    lapply(x$skeleton, `[[`, "is_mutable")
  )
  n <- x$n_chains
  x$skeleton <- lapply(seq_len(n), function(i) {
    i_int <- as.integer(i)
    list(
      times             = JuliaCall::julia_call("r_chain_times",                chains, chain = i_int),
      positions         = JuliaCall::julia_call("r_chain_positions",            chains, chain = i_int),
      velocities        = JuliaCall::julia_call("r_chain_velocities",           chains, chain = i_int),
      is_boomerang      = JuliaCall::julia_call("r_chain_is_boomerang",         chains, chain = i_int),
      is_mutable        = JuliaCall::julia_call("r_chain_is_mutable_boomerang", chains, chain = i_int),
      mu                = JuliaCall::julia_call("r_chain_mu",                   chains, chain = i_int)
    )
  })
  x$chains <- NULL
  x
}

#' @export
print.pdmp_result <- function(x, ...) {
  cat(sprintf(
    "PDMP result: %d dimension%s, %d chain%s\n",
    x$d, if (x$d == 1L) "" else "s",
    x$n_chains, if (x$n_chains == 1L) "" else "s"
  ))
  invisible(x)
}
