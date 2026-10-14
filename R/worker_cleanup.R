#' Release sampling state in a reused worker
#'
#' Clear package-owned Julia references left by sampling and request garbage
#' collection in both Julia and R. Use this after finishing a simulation task
#' in a persistent worker, before that worker starts another task. For example,
#' mirai daemons and foreach workers can process many models or data sets in
#' one R session; references held by JuliaCall can otherwise keep the previous
#' task's model and data alive after its R function returns.
#'
#' @details
#' Call this **inside the worker**, not in the parent R session: it only cleans
#' the process in which it runs. Register it with `on.exit(..., add = TRUE)`
#' at the start of the worker's task function so cleanup also runs after an
#' error. Finish extracting summaries or saving the fit before cleanup. Remove
#' large local fit objects when they are no longer needed; objects still
#' referenced by your code cannot be collected.
#'
#' This function does not start or shut down Julia, unload the package or
#' compiled libraries, delete streaming files, or remove saved results. It
#' keeps the worker ready for another task. Allocators can retain freed memory,
#' so a successful cleanup does not guarantee an immediate fall in the process's
#' reported memory use. Cleanup is usually unnecessary when each task's worker
#' exits immediately, or during ordinary interactive sampling where you still
#' want to use the fit.
#'
#' @md
#' @details
#' **Persistent worker example.**
#' The same task function can be submitted to mirai or foreach. Set the Julia
#' project in each worker when using a nondefault installation, and save a
#' materialized fit before returning a small result to the parent:
#'
#' ```r
#' run_condition <- function(seed, path, project_dir = NULL) {
#'   on.exit(PDMPSamplersR::cleanup_worker_state(), add = TRUE)
#'   PDMPSamplersR::pdmpsamplers_setup(verbose = FALSE,
#'                                   project_dir = project_dir)
#'   fit <- PDMPSamplersR::pdmp_sample(
#'     function(x) x, d = 2, T = 100, seed = seed,
#'     materialize = TRUE, show_progress = FALSE)
#'   saveRDS(fit, path)
#'   rm(fit)
#'   path
#' }
#'
#' # After starting mirai daemons:
#' job <- mirai::mirai(run_condition(seed, path),
#'   seed = 1L, path = tempfile(fileext = ".rds"),
#'   run_condition = run_condition)
#'
#' # Or, after registering a foreach process-worker backend:
#' # library(foreach)
#' # foreach(i = seq_len(10), .export = "run_condition") %dopar% {
#' #   run_condition(i, sprintf("fit-%03d.rds", i))
#' # }
#' ```
#'
#' Workers must be able to write the supplied paths. Use shared storage if the
#' parent needs to read files produced on another machine. Export or pass
#' `run_condition` and its inputs according to your scheduler's conventions.
#'
#' @return An invisible logical scalar: `TRUE` when the initialized bridge was
#'   cleaned, `FALSE` otherwise. Attributes `initialized`, `reason` and `error`
#'   distinguish `"cleaned"`, `"not_initialized"` (a harmless no-op), and
#'   `"error"` (cleanup failed). An uninitialized worker is not initialized as
#'   a side effect. Inspect these attributes when diagnosing worker cleanup;
#'   cleanup errors are reported in the result rather than replacing an error
#'   from the task being cleaned up.
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
