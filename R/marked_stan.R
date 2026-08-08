#' Compile and configure marked custom-Stan models
#'
#' `pdmp_subsample_hpp_path()` returns the installed thread-local subset-hook
#' header. `compile_pdmp_stan_model()` compiles a Stan source with the required
#' `--allow-undefined` and `USER_HEADER` settings.
#'
#' @param path Path to a Stan source file.
#' @return The header path, or the compiled shared-library path.
#' @export
pdmp_subsample_hpp_path <- function() {
  path <- system.file("stan", "pdmp_subsample.hpp", package = "PDMPSamplersR")
  if (!nzchar(path)) {
    candidate <- file.path("inst", "stan", "pdmp_subsample.hpp")
    if (file.exists(candidate)) path <- normalizePath(candidate)
  }
  if (!nzchar(path) || !file.exists(path)) {
    cli::cli_abort("The installed PDMP Stan subset header could not be found.")
  }
  path
}

#' @rdname pdmp_subsample_hpp_path
#' @export
compile_pdmp_stan_model <- function(path) {
  validate_type(path, type = "character", n = 1)
  if (!file.exists(path) || !grepl("\\.stan$", path)) {
    cli::cli_abort("Argument {.arg path} must identify an existing Stan source file.")
  }
  check_for_julia_setup()
  .pdmpsamplers_julia_call(
    "_compile_model_with_header",
    normalizePath(path), pdmp_subsample_hpp_path()
  )
}

#' Certified residual envelope for a marked Stan likelihood
#'
#' @param weights Nonnegative component-by-observation curvature weights.
#' @param growth_rates Optional nonnegative growth coefficient per component.
#' @return A declarative residual-envelope specification.
#' @export
stan_residual_envelope <- function(weights, growth_rates = NULL) {
  if (is.vector(weights)) weights <- matrix(as.numeric(weights), nrow = 1L)
  if (!is.matrix(weights) || !is.numeric(weights) || length(weights) == 0L ||
      any(!is.finite(weights)) || any(weights < 0)) {
    cli::cli_abort("Argument {.arg weights} must be a nonempty finite nonnegative numeric matrix.")
  }
  if (is.null(growth_rates)) growth_rates <- rep(0, nrow(weights))
  if (!is.numeric(growth_rates) || length(growth_rates) != nrow(weights) ||
      any(!is.finite(growth_rates)) || any(growth_rates < 0)) {
    cli::cli_abort("Argument {.arg growth_rates} must contain one finite nonnegative value per envelope component.")
  }
  structure(
    list(type = "separable", weights = unname(weights),
         growth_rates = as.numeric(growth_rates)),
    class = "stan_residual_envelope"
  )
}

#' Exact marked subsampling specification for a custom Stan model
#'
#' @param n_observations Total number of marked likelihood contributions.
#' @param subsample_size Number drawn without replacement at each proposal.
#' @param prior_standata Prior-only data list or JSON path for the same model.
#' @param residual_envelope A certified object from [stan_residual_envelope()]
#'   or a recognized package provider such as [omrf_residual_envelope()].
#' @param anchor Optional unconstrained anchor.
#' @param n_anchor_updates Number of warmup anchor updates (currently zero for
#'   generic custom-Stan providers).
#' @param use_anchor_bank Logical; reserved for analytic providers.
#' @param bank_capacity Positive anchor-bank capacity.
#' @return A validated marked-subsampling specification.
#' @export
stan_marked_subsampling <- function(n_observations, subsample_size,
                                    prior_standata, residual_envelope,
                                    anchor = NULL, n_anchor_updates = 0L,
                                    use_anchor_bank = FALSE,
                                    bank_capacity = 1L) {
  if (!rlang::is_integerish(n_observations, n = 1L, finite = TRUE) ||
      n_observations < 2L) {
    cli::cli_abort("Argument {.arg n_observations} must be an integer of at least 2.")
  }
  if (!rlang::is_integerish(subsample_size, n = 1L, finite = TRUE) ||
      subsample_size < 1L || subsample_size >= n_observations) {
    cli::cli_abort("Argument {.arg subsample_size} must satisfy 1 <= m < N.")
  }
  n_observations <- as.integer(n_observations)
  subsample_size <- as.integer(subsample_size)
  if (missing(prior_standata) || is.null(prior_standata)) {
    cli::cli_abort("Argument {.arg prior_standata} is required.")
  }
  if (!inherits(residual_envelope, "stan_residual_envelope")) {
    cli::cli_abort("Argument {.arg residual_envelope} must be a certified package residual-envelope object.")
  }
  if (ncol(residual_envelope$weights) != n_observations) {
    cli::cli_abort("Residual-envelope weights must have one column per observation.")
  }
  if (!rlang::is_integerish(n_anchor_updates, n = 1L, finite = TRUE) ||
      n_anchor_updates < 0L) {
    cli::cli_abort("Argument {.arg n_anchor_updates} must be nonnegative.")
  }
  n_anchor_updates <- as.integer(n_anchor_updates)
  if (n_anchor_updates > 0L) {
    cli::cli_abort("Generic custom-Stan marked sampling currently requires {.arg n_anchor_updates = 0L}.")
  }
  if (!is.logical(use_anchor_bank) || length(use_anchor_bank) != 1L || is.na(use_anchor_bank)) {
    cli::cli_abort("Argument {.arg use_anchor_bank} must be TRUE or FALSE.")
  }
  if (isTRUE(use_anchor_bank)) {
    cli::cli_abort("Generic custom-Stan marked sampling does not use an anchor bank.")
  }
  if (!rlang::is_integerish(bank_capacity, n = 1L, finite = TRUE) ||
      bank_capacity < 1L) {
    cli::cli_abort("Argument {.arg bank_capacity} must be positive.")
  }
  bank_capacity <- as.integer(bank_capacity)
  if (!is.null(anchor) && (!is.numeric(anchor) || any(!is.finite(anchor)))) {
    cli::cli_abort("Argument {.arg anchor} must be NULL or a finite numeric vector.")
  }
  structure(list(
    n_observations = n_observations, subsample_size = subsample_size,
    prior_standata = prior_standata, residual_envelope = residual_envelope,
    anchor = if (is.null(anchor)) NULL else as.numeric(anchor),
    n_anchor_updates = n_anchor_updates, use_anchor_bank = use_anchor_bank,
    bank_capacity = bank_capacity
  ), class = "stan_marked_subsampling")
}

#' OMRF person-level residual envelope
#'
#' Constructs the categorical-logit curvature bound
#' `0.5 * sum_j opnorm(B[n,j])^2` for each person. Threshold and interaction
#' arguments are Stan unconstrained block names and are resolved at sampler
#' initialization.
#'
#' @param X Integer person-by-node matrix encoded from zero.
#' @param seen Number of observed categories per node.
#' @param thresholds Threshold parameter block name.
#' @param interactions Interaction parameter block name; edge order is the
#'   upper triangle `(1,2), (1,3), ..., (P-1,P)`.
#' @param backend Selected-gradient backend. `"stan_header"` is currently
#'   implemented; `"analytic"` is reserved.
#' @return A certified Stan residual-envelope specification.
#' @export
omrf_residual_envelope <- function(X, seen, thresholds, interactions,
                                   backend = c("stan_header", "analytic")) {
  backend <- match.arg(backend)
  if (backend == "analytic") {
    cli::cli_abort("The analytic OMRF selected-gradient backend is not yet enabled; use {.val stan_header}.")
  }
  if (!is.matrix(X) || !is.numeric(X) || any(!is.finite(X)) || any(X != floor(X))) {
    cli::cli_abort("Argument {.arg X} must be a finite integer-valued person-by-node matrix.")
  }
  X <- matrix(as.integer(X), nrow(X), ncol(X))
  if (!rlang::is_integerish(seen, n = ncol(X), finite = TRUE) || any(seen < 1L)) {
    cli::cli_abort("Argument {.arg seen} must contain one positive category count per node.")
  }
  seen <- as.integer(seen)
  for (j in seq_len(ncol(X))) {
    if (any(X[, j] < 0L | X[, j] >= seen[j])) {
      cli::cli_abort("Column {j} of {.arg X} must be encoded in 0:(seen[j]-1).")
    }
  }
  if (!is.character(thresholds) || length(thresholds) != 1L ||
      !is.character(interactions) || length(interactions) != 1L) {
    cli::cli_abort("Arguments {.arg thresholds} and {.arg interactions} must be scalar block names.")
  }
  P <- ncol(X)
  threshold_starts <- cumsum(c(0L, seen[-length(seen)] - 1L))
  no_thresholds <- sum(seen - 1L)
  edges <- if (P > 1L) {
    do.call(rbind, lapply(seq_len(P - 1L), function(j) {
      cbind(j, seq.int(j + 1L, P))
    }))
  } else {
    matrix(integer(), nrow = 0L, ncol = 2L)
  }
  d_model <- no_thresholds + nrow(edges)
  weights <- numeric(nrow(X))
  for (n in seq_len(nrow(X))) {
    total <- 0
    for (j in seq_len(P)) {
      q <- seen[j] - 1L
      if (q == 0L) next
      B <- matrix(0, q, d_model)
      for (u in seq_len(q)) {
        B[u, threshold_starts[j] + u] <- 1
        incident <- which(edges[, 1L] == j | edges[, 2L] == j)
        neighbours <- ifelse(edges[incident, 1L] == j,
                             edges[incident, 2L], edges[incident, 1L])
        B[u, no_thresholds + incident] <- u * X[n, neighbours]
      }
      total <- total + norm(B, type = "2")^2
    }
    weights[n] <- 0.5 * total
  }
  envelope <- stan_residual_envelope(weights)
  envelope$type <- "omrf"
  envelope$X <- X
  envelope$seen <- seen
  envelope$thresholds <- thresholds
  envelope$interactions <- interactions
  envelope$edge_order <- edges
  class(envelope) <- c("omrf_residual_envelope", "stan_residual_envelope")
  envelope
}

.resolve_unconstrained_spec <- function(spec, unc_names, arg = "parameters") {
  if (is.null(spec)) return(seq_along(unc_names))
  if (rlang::is_integerish(spec)) {
    idx <- as.integer(spec)
    if (any(idx < 1L | idx > length(unc_names))) {
      cli::cli_abort("Integer {.arg {arg}} values must lie in 1:{length(unc_names)}.")
    }
  } else if (is.character(spec)) {
    idx <- unlist(lapply(spec, function(value) {
      exact <- which(unc_names == value)
      if (length(exact)) return(exact)
      prefix <- which(startsWith(unc_names, paste0(value, ".")) |
                      startsWith(unc_names, paste0(value, "[")))
      if (!length(prefix)) {
        cli::cli_abort("{.arg {arg}} value {.val {value}} did not match an unconstrained name or block prefix.")
      }
      prefix
    }), use.names = FALSE)
  } else {
    cli::cli_abort("Argument {.arg {arg}} must be NULL, character, or integer.")
  }
  if (anyDuplicated(idx)) cli::cli_abort("Resolved {.arg {arg}} coordinates contain duplicates.")
  idx
}

#' Inspect Stan unconstrained parameter mappings
#'
#' @param path_to_stanmodel Stan source or compiled-library path.
#' @param standata Named data list or JSON path.
#' @param parameters Optional unconstrained names, block prefixes, or indices.
#' @return A data frame with resolved 1-based indices and names.
#' @export
stan_parameter_mapping <- function(path_to_stanmodel, standata,
                                   parameters = NULL) {
  validate_type(path_to_stanmodel, type = "character", n = 1)
  data_path <- standata
  temporary <- FALSE
  if (is.list(standata)) {
    data_path <- tempfile(fileext = ".json")
    write_stan_json(standata, data_path)
    temporary <- TRUE
  }
  if (temporary) on.exit(unlink(data_path), add = TRUE)
  check_for_julia_setup()
  model_path <- normalizePath(path_to_stanmodel, mustWork = TRUE)
  if (grepl("\\.stan$", model_path)) {
    code <- paste(readLines(model_path, warn = FALSE), collapse = "\n")
    model_path <- if (grepl("pdmp_get_subsample_", code, fixed = TRUE)) {
      compile_pdmp_stan_model(model_path)
    } else {
      .pdmpsamplers_julia_call("_compile_model", model_path)
    }
  }
  unc_names <- .pdmpsamplers_julia_call(
    "r_get_param_unc_names", model_path, normalizePath(data_path, mustWork = TRUE)
  )
  idx <- .resolve_unconstrained_spec(parameters, unc_names)
  data.frame(index = idx, name = unc_names[idx], stringsAsFactors = FALSE)
}
