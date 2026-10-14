# Identity hashes used to check that a restored anchor bank belongs to the
# same data, parameter order and residual envelope (MD5 of the serialised
# object, or of the file contents for a data file).
.pdmp_object_hash <- function(value) {
  path <- tempfile("pdmpsamplers-flow-")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  writeBin(serialize(value, NULL, version = 3L), path)
  unname(tools::md5sum(path))
}

.pdmp_dataset_hash <- function(value) {
  if (is.character(value) && length(value) == 1L && file.exists(value)) {
    size <- file.info(value)$size
    con <- file(value, open = "rb")
    on.exit(close(con), add = TRUE)
    return(.pdmp_object_hash(readBin(con, what = "raw", n = size)))
  }
  .pdmp_object_hash(value)
}

#' Frozen position anchor bank from a set of states
#'
#' Builds the anchor bank used by position-selected OMRF subsampling from
#' states of a previous (typically full-gradient) chain, for example the online
#' grid snapshots of [read_online_summaries()]. It applies the same rule the
#' sampler uses for its own warmup: coordinate scales are the standard
#' deviations of the states (floored at a tenth of their median), the first
#' anchor is the last state, and further anchors are added by greedy
#' farthest-point selection in the scaled metric.
#'
#' During sampling each state uses the anchor nearest in a shared metric, so
#' the anchor is a function of position and the Voronoi cells are bounded by
#' hyperplanes. `metric = "scaled"` uses the coordinate scales above;
#' `metric = "hessian"` uses the mean likelihood Hessian of the anchors plus
#' the identity, which approximates the posterior precision and so predicts
#' the control-variate error better (it requires global HCV anchors).
#'
#' Pass the result as `subsampling_position_anchor_bank_spec` to
#' [pdmp_sample_from_stanmodel()], with `subsampling_anchor_selection_policy =
#' "position"`, `t_warmup = 0` and a frozen flow.
#'
#' @param standata The Stan data of the subsampled fit.
#' @param subsampling The [stan_subsampling()] object of the subsampled fit.
#' @param parameter_names Unconstrained parameter names, in sampler order.
#' @param states Numeric matrix with one state per column.
#' @param capacity Maximum number of anchors.
#' @param metric Anchor selection metric, `"scaled"` or `"hessian"`.
#' @return A version-1 position anchor bank specification.
#' @export
position_anchor_bank <- function(standata, subsampling, parameter_names,
                                 states, capacity = 8L,
                                 metric = c("scaled", "hessian")) {
  metric <- match.arg(metric)
  states <- as.matrix(states)
  if (nrow(states) != length(parameter_names) || ncol(states) < 1L) {
    cli::cli_abort("{.arg states} needs one row per parameter and at least one column.")
  }
  scales <- if (ncol(states) > 1L) apply(states, 1L, stats::sd) else
    rep(0, nrow(states))
  positive <- scales[scales > sqrt(.Machine$double.eps)]
  scales <- pmax(scales, if (length(positive)) 0.1 * stats::median(positive) else 1)

  k <- min(as.integer(capacity), ncol(states))
  chosen <- ncol(states)
  min_distance <- rep(Inf, ncol(states))
  while (TRUE) {
    last <- states[, chosen[length(chosen)]]
    min_distance <- pmin(min_distance, colSums(((states - last) / scales)^2))
    if (length(chosen) == k) break
    chosen <- c(chosen, which.max(min_distance))
  }

  list(version = 1L, selection_policy = "position",
       tie_policy = "lowest_index_outgoing", requested_capacity = k,
       ordered_anchor_positions = states[, chosen, drop = FALSE],
       distance_scales = scales, distance_metric = metric,
       parameter_order = as.character(parameter_names),
       dataset_identity = list(hash = .pdmp_dataset_hash(standata)),
       envelope_identity = list(hash = .pdmp_object_hash(subsampling$residual_envelope)))
}
