#' Clear condition-specific Julia/R state in a persistent worker.
#'
#' This is intentionally small and stable so external schedulers can release
#' BridgeStan and PDMP objects without knowing the JuliaCall Main roots used by
#' the sampling interface.
#'
#' @return An invisible logical scalar.  Its attributes distinguish a worker
#' that never initialized Julia from a worker whose initialized bridge could
#' not be cleaned.
#' @export
cleanup_worker_state <- function() {
  result <- function(value, initialized, reason, error = NULL) {
    structure(value, initialized = isTRUE(initialized), reason = reason,
              error = error)
  }
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    gc(full = TRUE)
    return(invisible(result(FALSE, FALSE, "not_initialized")))
  }
  tryCatch({
    julia_state <- get(".julia", envir = asNamespace("JuliaCall"), inherits = FALSE)
    if (!isTRUE(julia_state$initialized)) {
      result(FALSE, FALSE, "not_initialized")
    } else {
      bridge_loaded <- isTRUE(JuliaCall::julia_eval(
        "isdefined(Main, :PDMPSamplersRBridge)"))
      if (!bridge_loaded) {
        result(FALSE, TRUE, "error",
               "Julia is initialized but PDMPSamplersRBridge is unavailable.")
      } else {
        cleaned <- .pdmpsamplers_julia_call("r_cleanup_worker_state")
        if (isTRUE(cleaned) || is.null(cleaned)) {
          result(TRUE, TRUE, "cleaned")
        } else {
          result(FALSE, TRUE, "error",
                 "The Julia worker cleanup function returned a non-success value.")
        }
      }
    }
  }, error = function(e) {
    initialized <- tryCatch(isTRUE(get(".julia", envir =
      asNamespace("JuliaCall"), inherits = FALSE)$initialized),
      error = function(...) FALSE)
    result(FALSE, initialized, if (initialized) "error" else "not_initialized",
           conditionMessage(e))
  }) -> cleaned
  gc(full = TRUE)
  invisible(cleaned)
}
