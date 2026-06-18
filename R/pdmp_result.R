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
#' This only needs to be called manually when the result was created with
#' \code{materialize = FALSE}; all sampling functions materialise the trace
#' automatically by default.
#'
#' @param x A \code{pdmp_result}.
#' @returns A \code{pdmp_result} with \code{$chains = NULL} and
#'   \code{$skeleton} populated.
#' @export
materialize <- function(x) {
  .check_pdmp(x)
  chains <- if (!is.null(x$chains)) x$chains else .ensure_chains(x)
  n <- x$n_chains
  is_factorized <- .pdmpsamplers_julia_call("r_chain_is_factorized", chains, chain = 1L)
  x$skeleton <- if (is_factorized) {
    lapply(seq_len(n), function(i) {
      i_int <- as.integer(i)
      list(
        sparse             = TRUE,
        initial_time       = .pdmpsamplers_julia_call("r_chain_sparse_initial_time",      chains, chain = i_int),
        initial_position   = .pdmpsamplers_julia_call("r_chain_sparse_initial_position",  chains, chain = i_int),
        initial_velocity   = .pdmpsamplers_julia_call("r_chain_sparse_initial_velocity",  chains, chain = i_int),
        event_indices      = .pdmpsamplers_julia_call("r_chain_sparse_event_indices",     chains, chain = i_int),
        event_times        = .pdmpsamplers_julia_call("r_chain_sparse_event_times",       chains, chain = i_int),
        event_positions    = .pdmpsamplers_julia_call("r_chain_sparse_event_positions",   chains, chain = i_int),
        event_velocities   = .pdmpsamplers_julia_call("r_chain_sparse_event_velocities",  chains, chain = i_int),
        is_boomerang       = .pdmpsamplers_julia_call("r_chain_is_boomerang",             chains, chain = i_int),
        mu                 = .pdmpsamplers_julia_call("r_chain_mu",                       chains, chain = i_int)
      )
    })
  } else {
    lapply(seq_len(n), function(i) {
      i_int <- as.integer(i)
      list(
        sparse        = FALSE,
        times         = .pdmpsamplers_julia_call("r_chain_times",                chains, chain = i_int),
        positions     = .pdmpsamplers_julia_call("r_chain_positions",            chains, chain = i_int),
        velocities    = .pdmpsamplers_julia_call("r_chain_velocities",           chains, chain = i_int),
        is_boomerang  = .pdmpsamplers_julia_call("r_chain_is_boomerang",         chains, chain = i_int),
        is_mutable    = .pdmpsamplers_julia_call("r_chain_is_mutable_boomerang", chains, chain = i_int),
        mu            = .pdmpsamplers_julia_call("r_chain_mu",                   chains, chain = i_int)
      )
    })
  }
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

# Private alias used by sampling functions to call materialize() without the
# argument of the same name shadowing the function in their local scope.
.materialize <- function(x) materialize(x)
