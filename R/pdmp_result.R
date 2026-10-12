# Constructor (internal)
new_pdmp_result <- function(chains, stats, d, n_chains, skeleton = NULL,
                            initial_states = NULL, retained_initial_states = NULL,
                            terminal_states = NULL) {
  structure(
    list(
      chains   = chains,
      stats    = stats,
      d        = d,
      n_chains = n_chains,
      skeleton = skeleton,
      initial_states = initial_states,
      retained_initial_states = retained_initial_states,
      terminal_states = terminal_states
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
  if (is.null(x$chains) && !is.null(x$skeleton)) return(x)
  chains <- if (!is.null(x$chains)) x$chains else .ensure_chains(x)
  n <- x$n_chains
  has_live_states <- vapply(seq_len(n), function(i)
    isTRUE(.pdmpsamplers_julia_call("r_chain_has_live_state_provenance",
      chains, chain = as.integer(i))), logical(1L))
  x$initial_states <- if (all(has_live_states)) lapply(seq_len(n), function(i)
    .pdmpsamplers_julia_call("r_chain_initial_state",
      chains, chain = as.integer(i))) else NULL
  x$retained_initial_states <- if (all(has_live_states)) lapply(seq_len(n), function(i)
    .pdmpsamplers_julia_call("r_chain_retained_initial_state",
      chains, chain = as.integer(i))) else NULL
  x$terminal_states <- if (all(has_live_states)) lapply(seq_len(n), function(i)
    .pdmpsamplers_julia_call("r_chain_terminal_state",
      chains, chain = as.integer(i))) else NULL
  bank_specs <- attr(x, "position_anchor_bank_spec", exact = TRUE)
  if (!is.null(bank_specs) && !is.null(x$terminal_states)) {
    specs <- if (!is.null(bank_specs$version)) list(bank_specs) else bank_specs
    if (length(specs) == length(x$terminal_states)) {
      for (i in seq_along(specs)) {
        specs[[i]]$physical_endpoint <- x$terminal_states[[i]]
        specs[[i]]$frozen_flow <- x$terminal_states[[i]]$flow
      }
      attr(x, "position_anchor_bank_spec") <- if (length(specs) == 1L)
        specs[[1L]] else specs
    }
  }
  is_streaming <- .pdmpsamplers_julia_call(
    "r_chain_is_streaming", chains, chain = 1L)
  is_factorized <- .pdmpsamplers_julia_call(
    "r_chain_is_factorized", chains, chain = 1L)
  x$skeleton <- if (is_streaming) {
    lapply(seq_len(n), function(i) {
      i_int <- as.integer(i)
      list(
        streaming = TRUE,
        format = "pdmpsamplers_typed_stream_v1",
        directory = .pdmpsamplers_julia_call(
          "r_chain_streaming_directory", chains, chain = i_int),
        chunk_count = .pdmpsamplers_julia_call(
          "r_chain_streaming_chunk_count", chains, chain = i_int),
        output_bytes = .pdmpsamplers_julia_call(
          "r_chain_streaming_output_bytes", chains, chain = i_int),
        physical_events = .pdmpsamplers_julia_call(
          "r_chain_streaming_physical_events", chains, chain = i_int),
        computational_boundaries = .pdmpsamplers_julia_call(
          "r_chain_streaming_computational_boundaries", chains,
          chain = i_int),
        buffer_maximum = .pdmpsamplers_julia_call(
          "r_chain_streaming_buffer_maximum", chains, chain = i_int),
        is_boomerang = .pdmpsamplers_julia_call(
          "r_chain_is_boomerang", chains, chain = i_int),
        mu = .pdmpsamplers_julia_call(
          "r_chain_mu", chains, chain = i_int),
        flow_cov = .pdmpsamplers_julia_call(
          "r_chain_streaming_flow_cov", chains, chain = i_int)
      )
    })
  } else if (is_factorized) {
    lapply(seq_len(n), function(i) {
      i_int <- as.integer(i)
      list(
        sparse             = TRUE,
        initial_time       = .pdmpsamplers_julia_call("r_chain_sparse_initial_time",      chains, chain = i_int),
        initial_position   = .pdmpsamplers_julia_call("r_chain_sparse_initial_position",  chains, chain = i_int),
        initial_velocity   = .pdmpsamplers_julia_call("r_chain_sparse_initial_velocity",  chains, chain = i_int),
        initial_free_mask  = .pdmpsamplers_julia_call("r_chain_sparse_initial_free_mask", chains, chain = i_int),
        event_indices      = .pdmpsamplers_julia_call("r_chain_sparse_event_indices",     chains, chain = i_int),
        event_times        = .pdmpsamplers_julia_call("r_chain_sparse_event_times",       chains, chain = i_int),
        event_positions    = .pdmpsamplers_julia_call("r_chain_sparse_event_positions",   chains, chain = i_int),
        event_velocities   = .pdmpsamplers_julia_call("r_chain_sparse_event_velocities",  chains, chain = i_int),
        event_free_values  = .pdmpsamplers_julia_call("r_chain_sparse_event_free_values", chains, chain = i_int),
        is_boomerang       = .pdmpsamplers_julia_call("r_chain_is_boomerang",             chains, chain = i_int),
        mu                 = .pdmpsamplers_julia_call("r_chain_mu",                       chains, chain = i_int),
        flow_cov           = .pdmpsamplers_julia_call("r_chain_flow_cov",                 chains, chain = i_int)
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
        free_masks    = .pdmpsamplers_julia_call("r_chain_free_masks",           chains, chain = i_int),
        is_boomerang  = .pdmpsamplers_julia_call("r_chain_is_boomerang",         chains, chain = i_int),
        is_mutable    = .pdmpsamplers_julia_call("r_chain_is_mutable_boomerang", chains, chain = i_int),
        mu            = .pdmpsamplers_julia_call("r_chain_mu",                   chains, chain = i_int),
        flow_cov      = .pdmpsamplers_julia_call("r_chain_flow_cov",             chains, chain = i_int)
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
