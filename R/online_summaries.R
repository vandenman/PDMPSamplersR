#' Online summaries of the retained trajectory
#'
#' Configures exact continuous-time summaries that are accumulated while a
#' streaming chain runs. Instead of storing every event, the sampler keeps,
#' for `batches` equal physical-time batches of the retained phase, the
#' integrals of \eqn{x}, \eqn{x^2}, \eqn{x^4} and of the free (non-zero)
#' indicator for every coordinate, the integrals of the number of free
#' coordinates and its square, and the number of sticky transitions per
#' coordinate. Optionally it also records the state on a regular time grid.
#' Memory and disk use therefore do not grow with the number of events.
#'
#' Pass the result as `online_summaries` to [pdmp_sample()] or
#' [pdmp_sample_from_stanmodel()] together with `trace_storage = "streaming"`,
#' and read the summaries afterwards with [read_online_summaries()].
#'
#' @param batches Number of equal physical-time batches of the retained phase.
#' @param grid_spacing Physical time between state snapshots, or `NA` for none.
#' @param write_events If `FALSE` (default), no retained-phase events are
#'   written; the online summaries are the only retained output.
#' @param warmup_batches Number of equal physical-time batches of the warmup
#'   phase with the same integrals (no snapshots), or 0 (default) for none.
#'   They are read with `read_online_summaries(path, phase = "warmup")`.
#' @return A list of class `pdmp_online_summaries`.
#' @export
online_summary_options <- function(batches = 400L, grid_spacing = NA_real_,
                                   write_events = FALSE, warmup_batches = 0L) {
  if (!rlang::is_integerish(warmup_batches, n = 1L, finite = TRUE) ||
      warmup_batches < 0L) {
    cli::cli_abort("{.arg warmup_batches} must be a nonnegative integer.")
  }
  if (!rlang::is_integerish(batches, n = 1L, finite = TRUE) || batches < 1L) {
    cli::cli_abort("{.arg batches} must be a positive integer.")
  }
  if (length(grid_spacing) != 1L ||
      !(is.na(grid_spacing) || (is.finite(grid_spacing) && grid_spacing > 0))) {
    cli::cli_abort("{.arg grid_spacing} must be positive or {.code NA}.")
  }
  structure(list(batches = as.integer(batches),
                 grid_spacing = as.numeric(grid_spacing),
                 write_events = isTRUE(write_events),
                 warmup_batches = as.integer(warmup_batches)),
            class = "pdmp_online_summaries")
}

# The Julia bridge reads these settings from the environment when it builds
# the streaming trace storage; set them for the duration of one sampler call.
.set_online_summary_env <- function(online_summaries, trace_storage) {
  vars <- c("PDMPSAMPLERS_ONLINE_BATCHES", "PDMPSAMPLERS_ONLINE_GRID_SPACING",
            "PDMPSAMPLERS_STREAM_EVENTS", "PDMPSAMPLERS_ONLINE_WARMUP_BATCHES")
  old <- Sys.getenv(vars, unset = NA_character_, names = TRUE)
  restore <- function() {
    for (v in vars) {
      if (is.na(old[[v]])) Sys.unsetenv(v) else
        do.call(Sys.setenv, stats::setNames(list(old[[v]]), v))
    }
  }
  if (is.null(online_summaries)) {
    Sys.unsetenv(vars)
    return(restore)
  }
  if (!inherits(online_summaries, "pdmp_online_summaries")) {
    cli::cli_abort(
      "{.arg online_summaries} must come from {.fn online_summary_options}.")
  }
  if (!identical(trace_storage, "streaming")) {
    cli::cli_abort(
      "{.arg online_summaries} requires {.code trace_storage = \"streaming\"}.")
  }
  Sys.setenv(
    PDMPSAMPLERS_ONLINE_BATCHES = online_summaries$batches,
    PDMPSAMPLERS_ONLINE_GRID_SPACING = if (is.na(online_summaries$grid_spacing))
      "NaN" else format(online_summaries$grid_spacing, digits = 17),
    PDMPSAMPLERS_STREAM_EVENTS = if (online_summaries$write_events) "1" else "0",
    PDMPSAMPLERS_ONLINE_WARMUP_BATCHES = online_summaries$warmup_batches %||% 0L)
  restore
}

#' Read online summaries written by a streaming chain
#'
#' @param path Either the `online_summaries.bin` file or a streaming trace
#'   directory, which is searched for exactly one such file of `phase`.
#' @param phase `"main"` (default) for the retained phase or `"warmup"`;
#'   only used when `path` is a directory.
#' @return A list with the batch integrals (`x`, `x2`, `x4`, `free`: one column
#'   per batch), `batch_time`, `size`, `size2` (integrals of the number of free
#'   coordinates and its square), `transitions` (per coordinate and batch),
#'   the snapshot grid (`grid_times`, `grid_x`, `grid_free`), and
#'   `start_time`, `width` and `end_time` of the retained phase.
#' @export
read_online_summaries <- function(path, phase = c("main", "warmup")) {
  phase <- match.arg(phase)
  if (dir.exists(path)) {
    found <- list.files(path, "^online_summaries\\.bin$", recursive = TRUE,
                        full.names = TRUE)
    found <- found[basename(dirname(found)) == phase]
    if (length(found) != 1L) {
      cli::cli_abort("Expected one online_summaries.bin under {.path {path}}, found {length(found)}.")
    }
    path <- found
  }
  con <- file(path, "rb")
  on.exit(close(con))
  if (!identical(rawToChar(readBin(con, "raw", 8L)), "PDMPONL1")) {
    cli::cli_abort("{.path {path}} is not an online summary file.")
  }
  dims <- readBin(con, "integer", 3L, size = 8L)
  d <- dims[[1L]]; nb <- dims[[2L]]; ng <- dims[[3L]]
  times <- readBin(con, "double", 3L)
  dbl <- function(n) readBin(con, "double", n)
  per_batch <- function() matrix(dbl(d * nb), d, nb)
  out <- list(d = d, n_batches = nb, start_time = times[[1L]],
              width = times[[2L]], end_time = times[[3L]])
  out$batch_time <- dbl(nb)
  out$x <- per_batch()
  out$x2 <- per_batch()
  out$x4 <- per_batch()
  out$free <- per_batch()
  out$size <- dbl(nb)
  out$size2 <- dbl(nb)
  out$transitions <- matrix(readBin(con, "integer", d * nb, size = 8L), d, nb)
  out$grid_times <- dbl(ng)
  out$grid_x <- matrix(dbl(d * ng), d, ng)
  out$grid_free <- matrix(as.logical(readBin(con, "raw", d * ng)), d, ng)
  out
}
