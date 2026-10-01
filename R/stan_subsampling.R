#' Compile and configure subsampling custom-Stan models
#'
#' `pdmp_subsample_hpp_path()` returns the installed thread-local subset-hook
#' header. `compile_pdmp_stan_model()` compiles a Stan source with the required
#' `--allow-undefined` and `USER_HEADER` settings.
#'
#' @param path Path to a Stan source file.
#' @return The header path, or the compiled shared-library path.
#' @name stan_subsampling
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

#' @rdname stan_subsampling
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

#' Certified residual envelope for a subsampling Stan likelihood
#'
#' @param weights Nonnegative component-by-observation curvature weights.
#' @param growth_rates Optional nonnegative growth coefficient per component.
#' @return A declarative residual-envelope specification.
#' @rdname stan_subsampling
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

#' Exact observation-subsampling specification for a custom Stan model
#'
#' @param n_observations Total number of subsampling likelihood contributions.
#' @param subsample_size Number drawn without replacement at each proposal.
#' @param prior_standata Prior-only data list or JSON path for the same model.
#' @param residual_envelope A certified object from [stan_residual_envelope()]
#'   or a recognized package provider such as [omrf_residual_envelope()].
#' @param anchor Optional unconstrained anchor.
#' @return A validated observation-subsampling specification.
#' @rdname stan_subsampling
#' @export
stan_subsampling <- function(n_observations, subsample_size,
                                    prior_standata, residual_envelope,
                                    anchor = NULL) {
  if (!rlang::is_integerish(n_observations, n = 1L, finite = TRUE) ||
      n_observations < 2L) {
    cli::cli_abort("Argument {.arg n_observations} must be an integer of at least 2.")
  }
  if (!rlang::is_integerish(subsample_size, n = 1L, finite = TRUE) ||
      subsample_size < 1L || subsample_size > n_observations) {
    cli::cli_abort("Argument {.arg subsample_size} must satisfy 1 <= m <= N.")
  }
  n_observations <- as.integer(n_observations)
  subsample_size <- as.integer(subsample_size)
  if (missing(prior_standata) || is.null(prior_standata)) {
    cli::cli_abort("Argument {.arg prior_standata} is required.")
  }
  if (!inherits(residual_envelope, "stan_residual_envelope")) {
    cli::cli_abort("Argument {.arg residual_envelope} must be a certified package residual-envelope object.")
  }
  if (inherits(residual_envelope, "omrf_residual_envelope")) {
    n_persons <- residual_envelope$n_persons
    n_nodes <- residual_envelope$n_nodes
    factorization <- residual_envelope$factorization
    if (!rlang::is_integerish(n_persons, n = 1L, finite = TRUE) ||
        n_persons < 1L ||
        !rlang::is_integerish(n_nodes, n = 1L, finite = TRUE) ||
        n_nodes < 1L ||
        !is.character(factorization) || length(factorization) != 1L ||
        !factorization %in% c("person", "person_node")) {
      cli::cli_abort(paste0(
        "Specialized OMRF envelopes require positive integer {.field n_persons} ",
        "and {.field n_nodes} metadata and a supported {.field factorization}."))
    }
    expected <- if (factorization == "person") {
      as.integer(n_persons)
    } else {
      as.integer(n_persons) * as.integer(n_nodes)
    }
    if (n_observations != expected) {
      cli::cli_abort(paste0(
        "OMRF {.field ", factorization, "} factorization requires ", expected,
        " likelihood contributions; received ", n_observations, "."))
    }
    compact_norm <- factorization == "person_node" && isTRUE(
      residual_envelope$bound_type %in%
        c("node_local_norm_low_rank", "norm_local_hcv_compact"))
    if (compact_norm) {
      powers <- if (identical(residual_envelope$bound_type,
                              "norm_local_hcv_compact")) 4L else 3L
      weights <- residual_envelope$norm_weights
      if (!is.matrix(weights) || !is.numeric(weights) ||
          !identical(dim(weights), c(powers * as.integer(n_nodes),
                                     as.integer(n_persons))) ||
          any(!is.finite(weights)) || any(weights < 0)) {
        cli::cli_abort(paste0(
          "Compact OMRF norm weights must be a finite nonnegative ",
          powers, "P by N matrix matching {.field n_nodes} and ",
          "{.field n_persons}."))
      }
    }
    if (factorization == "person_node" &&
        residual_envelope$bound_type %in%
          c("node_box_first_order", "cell_marked_node_box_first_order")) {
      membership <- residual_envelope$node_box_membership
      nodes <- residual_envelope$node_box_nodes
      minima <- residual_envelope$node_box_minima
      maxima <- residual_envelope$node_box_maxima
      if (!is.matrix(membership) ||
          !identical(dim(membership), c(as.integer(n_nodes),
                                        as.integer(n_persons))) ||
          !is.numeric(membership) || any(!is.finite(membership)) ||
          any(membership != floor(membership)) || !length(nodes) ||
          any(membership < 1L | membership > length(nodes)) ||
          !is.matrix(minima) || !is.matrix(maxima) ||
          !identical(dim(minima), c(length(nodes), as.integer(n_nodes))) ||
          !identical(dim(maxima), dim(minima)) ||
          any(!is.finite(minima)) || any(!is.finite(maxima)) ||
          any(minima > maxima)) {
        cli::cli_abort(paste0(
          "Compact OMRF node-box specifications require P-by-N integer ",
          "membership and matching finite component-by-P extrema."))
      }
    }
  } else if (is.null(residual_envelope$weights) ||
             ncol(residual_envelope$weights) != n_observations) {
    cli::cli_abort("Residual-envelope weights must have one column per observation.")
  }
  if (!is.null(anchor) && (!is.numeric(anchor) || any(!is.finite(anchor)))) {
    cli::cli_abort("Argument {.arg anchor} must be NULL or a finite numeric vector.")
  }
  if (inherits(residual_envelope, "omrf_residual_envelope") &&
      identical(residual_envelope$prior_backend, "analytic")) {
    residual_envelope <- .attach_omrf_analytic_prior(
      residual_envelope, prior_standata, "prior_standata")
  }
  structure(list(
    n_observations = n_observations, subsample_size = subsample_size,
    prior_standata = prior_standata, residual_envelope = residual_envelope,
    anchor = if (is.null(anchor)) NULL else as.numeric(anchor)
  ), class = "stan_subsampling")
}

.attach_omrf_analytic_prior <- function(spec, standata, argument = "standata") {
  if (!is.list(standata)) {
    cli::cli_abort("Analytic OMRF priors require {.arg {argument}} as a named list.")
  }
  required <- c("prior_threshold_alpha", "prior_threshold_beta")
  if (!all(required %in% names(standata))) {
    cli::cli_abort("Analytic OMRF priors require threshold alpha and beta data.")
  }
  # The OMRF simulation data can intentionally carry both independent-prior
  # and node-scale fields so that the same generated data object can be routed
  # to either Stan model.  For the native full-gradient independent model,
  # `prior_cauchy_scale` is the model-defining field and must take precedence;
  # node-scale prior-only data are stripped by `make_prior_stan_data()` before
  # this helper is called for the dependent route.
  interaction <- if (!is.null(standata$prior_cauchy_scale)) {
    list(type = "cauchy", scale = standata$prior_cauchy_scale)
  } else if (!is.null(standata$slab_base_scale) &&
             !is.null(standata$prior_global_scale_sd) &&
             !is.null(standata$prior_node_scale_sd)) {
    list(
      type = "node_scale_gaussian",
      scale = standata$slab_base_scale,
      global_scale_sd = standata$prior_global_scale_sd,
      node_scale_sd = standata$prior_node_scale_sd
    )
  } else if (!is.null(standata$prior_interaction_sd)) {
    list(type = "gaussian", scale = standata$prior_interaction_sd)
  } else {
    cli::cli_abort(paste0(
      "The requested analytic OMRF prior is not a supported independent ",
      "Cauchy/Gaussian or node-scale Gaussian prior."))
  }
  values <- c(standata[required], interaction$scale,
              interaction$global_scale_sd, interaction$node_scale_sd)
  if (any(!vapply(values, function(value) {
    is.numeric(value) && length(value) == 1L && is.finite(value) && value > 0
  }, logical(1)))) {
    cli::cli_abort("Analytic OMRF prior parameters must be finite and positive.")
  }
  spec$analytic_prior <- list(
    type = interaction$type,
    threshold_alpha = standata$prior_threshold_alpha,
    threshold_beta = standata$prior_threshold_beta,
    interaction_scale = interaction$scale,
    global_scale = if (identical(interaction$type, "node_scale_gaussian"))
      "global_scale" else NULL,
    node_scale = if (identical(interaction$type, "node_scale_gaussian"))
      "node_scale" else NULL,
    global_scale_sd = interaction$global_scale_sd,
    node_scale_sd = interaction$node_scale_sd)
  spec
}

#' OMRF person-level residual envelope
#'
#' Constructs a stacked categorical-logit curvature bound for each person.
#' If `B_n` vertically stacks the node-specific affine logit maps, the block
#' diagonal categorical covariance has operator norm at most `1/2`, hence the
#' complete person Hessian is bounded by
#' `0.5 * opnorm(B_n)^2`. This is no larger than summing the separate node
#' spectral bounds and preserves a global, anchor-independent certificate.
#' The analytic backend uses a tighter low-rank node-local certificate. For
#' each node it expands the product of threshold and interaction terms against
#' the Euclidean norm of that person's neighbour-category vector. This needs
#' only three envelope components per node while retaining covariate magnitude
#' and preventing thresholds or nonincident interactions from inflating the
#' bound.
#' Threshold and interaction arguments are Stan unconstrained block names and
#' are resolved at sampler initialization.
#'
#' @param X Integer person-by-node matrix encoded from zero.
#' @param seen Number of observed categories per node.
#' @param thresholds Threshold parameter block name.
#' @param interactions Interaction parameter block name; edge order is the
#'   upper triangle `(1,2), (1,3), ..., (P-1,P)`.
#' @param backend Residual-gradient backend. The default `"analytic"` evaluates
#'   the exact OMRF likelihood residual without invoking Stan autodiff.
#'   `"stan"` retains the selected-likelihood Stan path for validation.
#' @param use_hcv Whether to use the analytic categorical Hessian control
#'   variate. This requires `backend = "analytic"` and retains a certified
#'   first-order fallback through damping.
#' @param hcv_damping Positive damping scale. By default this is the number of
#'   OMRF likelihood coordinates. `Inf` selects the exact undamped Hessian
#'   control variate and is currently supported only with
#'   `geometry = "pattern_local"`, `geometry = "grouped_box"`, or one of
#'   the node-grouped geometries.
#' @param hcv_scope HCV damping scope. The default `"global"` uses one
#'   displacement scale for the complete OMRF likelihood. `"node"` uses a
#'   separate displacement and Hessian action for each OMRF node, so a distant
#'   interaction only damps the HCV correction for its incident nodes.
#' @param zigzag_sign_aware Whether to use coordinatewise deterministic-sign
#'   subset bounds for Zig-Zag. This is an experimental analytic-OMRF option.
#' @param hcv_warmup HCV policy during warmup. `"hcv"` preserves the existing
#'   behavior. Experimental `"first_order"` uses the certified first-order
#'   residual bound during warmup and activates HCV at the main-phase boundary.
#' @param hcv_bank_learning Research-only frozen-bank learning policy.
#'   `"geometric"` retains standardized position coverage. `"exact_variance"`
#'   uses time-uniform states from each chain's own full-gradient warmup and
#'   greedily minimizes the exact finite-population HCV residual variance under
#'   the final adapted velocity metric. All learning and preparation are
#'   charged to that chain's warmup; the bank is frozen before retained sampling.
#'   `"selector_coupled_variance"` is a research-only variant which minimizes
#'   that same variance after assigning every warmup state with the actual
#'   standardized-position selector used during retained sampling.
#' @param factorization Experimental likelihood contribution layout. `"person"`
#'   retains one contribution per person. `"person_node"` exposes each
#'   person-node conditional as a separate sparse contribution.
#' @param factor_sampling Factor-subset design. `"size_biased"` uses the
#'   generic contribution design; `"node_stratified"` draws an equal number
#'   from every node while retaining exact size-biased proposal mixing.
#' @param direct_cv_cell_scope Research-only scope for
#'   `geometry = "direct_person_node_cv"`. `"cell"` constructs one residual
#'   table per deterministic grid cell. `"event_horizon"` constructs one
#'   certified table over the complete finite event-search horizon and reuses
#'   it for every internal deterministic cell.
#' @param prior_backend Deterministic prior-gradient backend. `"stan"` supports
#'   arbitrary custom priors. `"analytic"` recognizes the package OMRF
#'   independent Gaussian and Cauchy prior contracts from `prior_standata`.
#' @param geometry Certified residual-rate geometry. `"norm_local"` uses a
#'   compact covariate-norm decomposition. `"covariate_local"` retains
#'   empirical neighbour covariates in threshold, edge, and edge-pair terms
#'   while sharing trajectory scales across people. `"pattern_local"` groups
#'   identical node-neighbour patterns and retains their signed directional
#'   logit geometry, which is especially tighter for scalar BPS/Boomerang
#'   rates but is more expensive to construct.
#'   `"grouped_box"` is the preferred exact HCV clock for non-factorized
#'   subsampling: setup partitions whole persons into a fixed number of
#'   covariate boxes and retained sampling evaluates one certificate per box,
#'   never one certificate per person. It currently requires
#'   `backend = "analytic"`, `use_hcv = TRUE`, `hcv_damping = Inf`, and
#'   `factorization = "person"`.
#'   `"node_grouped_box"` independently partitions each node by its incident
#'   covariates and reuses the node's trajectory geometry across those boxes.
#'   Experimental `"cell_cached_node_grouped_box"` uses the same node boxes
#'   but samples marked persons from the already-certified grid-cell masses,
#'   then thins only the selected person using its exact observed covariate
#'   pattern. This avoids a second all-box pointwise traversal at every clock
#'   proposal while retaining a whole-person stochastic gradient.
#'   Experimental `"lazy_node_grouped_box"` builds the same recursive leaves
#'   but exposes a coarser root frontier to GridThinning. At an accepted point
#'   proposal it descends the certified support tree by exact thinning, so
#'   child and leaf certificates are evaluated only along the sampled branch.
#' @param person_signed_roof Research-only linear BPS certificate policy for
#'   `geometry = "direct_factor_roof_person_hcv"`. Signed node residual
#'   endpoints are summed within each person before taking an absolute bound;
#'   the old factor-summed certificate remains a valid fallback. The estimator
#'   and retained target are unchanged. Other flows keep the existing roof.
#' @param signed_only_person_roof Research-only linear BPS policy which builds
#'   the certified signed whole-person endpoint roof without also building
#'   the factor-summed fallback table. It preserves the estimator and uses
#'   pointwise factor refinement for the selected person.
#' @param joint_signed_group_size Research-only linear BPS joint signed
#'   marked roof. Use `4` for fixed seeded groups of four, `16` for the
#'   opt-in N=500, P=15, K=4 weighted ragged-group research route, or `0`
#'   to keep the ordinary scheduler. Requires `signed_only_person_roof = TRUE`.
#' @param shared_current_categorical Opt-in JS4/JS16 evaluator that reuses
#'   current categorical moments for exact active-neighbour patterns. Supports
#'   the independent-prior linear BPS research route only.
#' @param box_root_groups Root boxes per node used by
#'   `geometry = "lazy_node_grouped_box"`. This must not exceed `box_groups`.
#' @param box_groups Maximum number of deterministic setup-time covariate
#'   boxes used by `geometry = "grouped_box"`, or boxes per node for
#'   `geometry = "node_grouped_box"`. The actual number is capped by the
#'   number of distinct relevant person patterns.
#' @return A certified Stan residual-envelope specification.
#' @rdname stan_subsampling
#' @export
.omrf_covariate_blocks <- function(X, incidence) {
  P <- ncol(X)
  N <- nrow(X)
  blocks <- vector("list", P)
  for (node in seq_len(P)) {
    neighbours <- incidence[[node]]$neighbours
    degree <- length(neighbours)
    node_rows <- matrix(1, nrow = 1L + 2L * degree + degree^2,
                        ncol = N)
    if (degree) {
      node_covariates <- X[, neighbours, drop = FALSE]
      node_rows[1L + seq_len(degree), ] <- t(node_covariates)
      node_rows[1L + degree + seq_len(degree), ] <- t(node_covariates)
      destination <- 1L + 2L * degree
      for (left_edge in seq_len(degree)) {
        for (right_edge in seq_len(degree)) {
          destination <- destination + 1L
          node_rows[destination, ] <- node_covariates[, left_edge] *
            node_covariates[, right_edge]
        }
      }
    }
    blocks[[node]] <- node_rows
  }
  blocks
}

.omrf_node_support_boxes <- function(X, incidence, group_limit) {
  split_boxes <- function(features) {
    if (!ncol(features)) {
      return(list(membership = rep.int(1L, nrow(features)), boxes = list(
        seq_len(nrow(features)))))
    }
    boxes <- list(seq_len(nrow(features)))
    while (length(boxes) < group_limit) {
      candidates <- lapply(boxes, function(rows) {
        ranges <- apply(features[rows, , drop = FALSE], 2L, function(x) {
          diff(range(x))
        })
        list(score = sum(ranges) * length(rows), ranges = ranges)
      })
      splittable <- which(vapply(
        candidates, function(x) any(x$ranges > 0), logical(1L)))
      if (!length(splittable)) break
      selected <- splittable[[which.max(vapply(
        candidates[splittable], `[[`, numeric(1L), "score"))]]
      rows <- boxes[[selected]]
      dimension <- which.max(candidates[[selected]]$ranges)
      ordered <- rows[order(features[rows, dimension], rows)]
      changes <- which(diff(features[ordered, dimension]) > 0)
      cut <- changes[[which.min(abs(changes - length(rows) / 2))]]
      boxes[[selected]] <- ordered[seq_len(cut)]
      boxes[[length(boxes) + 1L]] <- ordered[seq.int(cut + 1L,
                                                     length(ordered))]
    }
    membership <- integer(nrow(features))
    for (group in seq_along(boxes)) membership[boxes[[group]]] <- group
    list(membership = membership, boxes = boxes)
  }

  P <- ncol(X)
  N <- nrow(X)
  membership <- matrix(0L, nrow = P, ncol = N)
  nodes <- integer()
  minima <- list()
  maxima <- list()
  offset <- 0L
  for (node in seq_len(P)) {
    neighbours <- incidence[[node]]$neighbours
    features <- X[, neighbours, drop = FALSE]
    partition <- split_boxes(features)
    membership[node, ] <- offset + partition$membership
    for (group in seq_along(partition$boxes)) {
      lower <- upper <- numeric(P)
      if (length(neighbours)) {
        rows <- partition$boxes[[group]]
        lower[neighbours] <- apply(features[rows, , drop = FALSE], 2L, min)
        upper[neighbours] <- apply(features[rows, , drop = FALSE], 2L, max)
      }
      minima[[offset + group]] <- lower
      maxima[[offset + group]] <- upper
    }
    nodes <- c(nodes, rep.int(node, length(partition$boxes)))
    offset <- offset + length(partition$boxes)
  }
  list(membership = membership, nodes = nodes,
       minima = do.call(rbind, minima), maxima = do.call(rbind, maxima))
}

omrf_residual_envelope <- function(X, seen, thresholds, interactions,
                                   backend = c("analytic", "stan"),
                                   use_hcv = FALSE, hcv_damping = NULL,
                                   hcv_scope = c("global", "node"),
                                   zigzag_sign_aware = FALSE,
                                   hcv_warmup = c("hcv", "first_order"),
                                   hcv_bank_learning = c("geometric", "exact_variance",
                                                         "selector_coupled_variance"),
                                   factorization = c("person", "person_node"),
                                   factor_sampling = c("size_biased", "node_stratified"),
                                   direct_cv_cell_scope = c("cell", "event_horizon"),
                                   prior_backend = c("stan", "analytic"),
                                   exact_selected_refinement = FALSE,
                                   categorical_residual_cap = FALSE,
                                   signed_bps = FALSE,
                                   signed_harmonic = FALSE,
                                   person_signed_roof = FALSE,
                                   signed_only_person_roof = FALSE,
                                   joint_signed_group_size = 0L,
                                   shared_current_categorical = FALSE,
                                   geometry = c("norm_local", "norm_local_hcv_compact",
                                                "node_box_first_order",
                                                "cell_marked_node_box_first_order",
                                                "direct_person_node_hcv",
                                                "direct_factor_roof_person_hcv",
                                                "direct_person_node_cv",
                                                "covariate_local",
                                                "covariate_local_compact",
                                                "pattern_local", "grouped_box",
                                                "node_grouped_box",
                                                "cell_cached_node_grouped_box",
                                                "lazy_node_grouped_box"),
                                   box_groups = 16L,
                                   box_root_groups = max(1L, box_groups %/% 2L)) {
  backend <- match.arg(backend)
  hcv_scope <- match.arg(hcv_scope)
  hcv_warmup <- match.arg(hcv_warmup)
  hcv_bank_learning <- match.arg(hcv_bank_learning)
  factorization <- match.arg(factorization)
  factor_sampling <- match.arg(factor_sampling)
  direct_cv_cell_scope <- match.arg(direct_cv_cell_scope)
  prior_backend <- match.arg(prior_backend)
  geometry <- match.arg(geometry)
  if (!rlang::is_integerish(box_groups, n = 1L, finite = TRUE) ||
      box_groups < 1L) {
    cli::cli_abort("Argument {.arg box_groups} must be one positive integer.")
  }
  box_groups <- as.integer(box_groups)
  if (!rlang::is_integerish(box_root_groups, n = 1L, finite = TRUE) ||
      box_root_groups < 1L || box_root_groups > box_groups) {
    cli::cli_abort(paste0(
      "Argument {.arg box_root_groups} must be one positive integer no ",
      "larger than {.arg box_groups}."))
  }
  box_root_groups <- as.integer(box_root_groups)
  if (!is.logical(use_hcv) || length(use_hcv) != 1L || is.na(use_hcv)) {
    cli::cli_abort("Argument {.arg use_hcv} must be TRUE or FALSE.")
  }
  if (!is.logical(zigzag_sign_aware) || length(zigzag_sign_aware) != 1L ||
      is.na(zigzag_sign_aware)) {
    cli::cli_abort("Argument {.arg zigzag_sign_aware} must be TRUE or FALSE.")
  }
  if (!is.logical(exact_selected_refinement) ||
      length(exact_selected_refinement) != 1L ||
      is.na(exact_selected_refinement)) {
    cli::cli_abort(
      "Argument {.arg exact_selected_refinement} must be TRUE or FALSE.")
  }
  if (!is.logical(categorical_residual_cap) ||
      length(categorical_residual_cap) != 1L ||
      is.na(categorical_residual_cap)) {
    cli::cli_abort(
      "Argument {.arg categorical_residual_cap} must be TRUE or FALSE.")
  }
  if (!is.logical(signed_bps) || length(signed_bps) != 1L ||
      is.na(signed_bps)) {
    cli::cli_abort("Argument {.arg signed_bps} must be TRUE or FALSE.")
  }
  if (!is.logical(signed_harmonic) || length(signed_harmonic) != 1L ||
      is.na(signed_harmonic)) {
    cli::cli_abort("Argument {.arg signed_harmonic} must be TRUE or FALSE.")
  }
  if (!is.logical(person_signed_roof) || length(person_signed_roof) != 1L ||
      is.na(person_signed_roof)) {
    cli::cli_abort("Argument {.arg person_signed_roof} must be TRUE or FALSE.")
  }
  if (person_signed_roof && geometry != "direct_factor_roof_person_hcv") {
    cli::cli_abort(paste0(
      "Research {.arg person_signed_roof} requires ",
      "{.arg geometry = \"direct_factor_roof_person_hcv\"}."))
  }
  if (!is.logical(signed_only_person_roof) ||
      length(signed_only_person_roof) != 1L ||
      is.na(signed_only_person_roof)) {
    cli::cli_abort("Argument {.arg signed_only_person_roof} must be TRUE or FALSE.")
  }
  if (signed_only_person_roof &&
      !(person_signed_roof && geometry == "direct_factor_roof_person_hcv" &&
        direct_cv_cell_scope == "cell")) {
    cli::cli_abort(paste0(
      "Research {.arg signed_only_person_roof} requires signed whole-person ",
      "HCV roofs and per-cell certificates."))
  }
  if (!rlang::is_integerish(joint_signed_group_size, n = 1L,
                           finite = TRUE) ||
      !joint_signed_group_size %in% c(0L, 4L, 16L, 32L, 64L)) {
    cli::cli_abort("Research {.arg joint_signed_group_size} must be 0, 4, 16, 32, or 64.")
  }
  if (joint_signed_group_size == 4L &&
      (!signed_only_person_roof || nrow(X) %% 4L != 0L)) {
    cli::cli_abort(paste0("Research joint groups of four require a ",
      "divisible signed-only person-HCV population."))
  }
  if (!is.logical(shared_current_categorical) ||
      length(shared_current_categorical) != 1L ||
      is.na(shared_current_categorical)) {
    cli::cli_abort("Argument {.arg shared_current_categorical} must be TRUE or FALSE.")
  }
  if (shared_current_categorical &&
      joint_signed_group_size == 0L) {
    cli::cli_abort("Shared-current categorical evaluation requires joint groups.")
  }
  if (joint_signed_group_size >= 16L &&
      !(shared_current_categorical && signed_only_person_roof &&
        backend == "analytic" && isTRUE(use_hcv) &&
        is.infinite(hcv_damping) && hcv_scope == "global" &&
        geometry == "direct_factor_roof_person_hcv" &&
        factorization == "person" && prior_backend == "analytic" &&
        nrow(X) >= 2L * joint_signed_group_size && all(seen == 4L))) {
    cli::cli_abort(paste0("Grouped joint signed bounds require N >= 2m, four categories per node, ",
      "shared-current signed global undamped HCV with analytic prior."))
  }
  if (isTRUE(signed_bps) &&
      (backend != "analytic" || factorization != "person_node" ||
       isTRUE(use_hcv) || geometry != "covariate_local_compact")) {
    cli::cli_abort(paste0(
      "Research {.arg signed_bps} requires analytic, no-HCV, person-node, ",
      "covariate_local_compact OMRF controls."))
  }
  if (isTRUE(use_hcv) && backend != "analytic") {
    cli::cli_abort("OMRF Hessian control variates require {.arg backend = \"analytic\"}.")
  }
  if (hcv_scope == "node" && !isTRUE(use_hcv)) {
    cli::cli_abort("Node-local HCV requires {.arg use_hcv = TRUE}.")
  }
  if (factorization == "person_node" && backend != "analytic") {
    cli::cli_abort("Person-node OMRF factorization requires {.arg backend = \"analytic\"}.")
  }
  if (factorization == "person" && factor_sampling != "size_biased") {
    cli::cli_abort("Node-stratified sampling requires {.arg factorization = \"person_node\"}.")
  }
  if (geometry == "covariate_local_compact" &&
      factorization != "person_node") {
    cli::cli_abort("Compact covariate-local envelopes require {.arg factorization = \"person_node\"}.")
  }
  if (geometry %in% c("norm_local_hcv_compact", "direct_person_node_hcv") &&
      !(backend == "analytic" && isTRUE(use_hcv) &&
        is.infinite(hcv_damping) && hcv_scope == "global" &&
        hcv_warmup == "hcv" &&
        factorization == "person_node" &&
        factor_sampling == "size_biased")) {
    cli::cli_abort(paste0(
      "Research HCV geometry requires ",
      "analytic, undamped global HCV, person-node factorization, and ",
      "size-biased factor sampling, active from warmup."))
  }
  if (geometry == "direct_factor_roof_person_hcv" &&
      !(backend == "analytic" && isTRUE(use_hcv) &&
        is.infinite(hcv_damping) && hcv_scope == "global" &&
        hcv_warmup == "hcv" && factorization == "person" &&
        factor_sampling == "size_biased")) {
    cli::cli_abort(paste0(
      "Research factor-roof person HCV geometry requires analytic, ",
      "undamped global HCV, person factorization, size-biased person ",
      "sampling, active from warmup."))
  }
  if (geometry == "direct_person_node_cv" &&
      !(backend == "analytic" && !isTRUE(use_hcv) &&
        factorization == "person_node" &&
        factor_sampling == "size_biased")) {
    cli::cli_abort(paste0(
      "Research direct first-order geometry requires analytic, no-HCV, ",
      "person-node factorization and size-biased factor sampling."))
  }
  if (geometry %in% c("node_box_first_order",
                      "cell_marked_node_box_first_order") &&
      !(backend == "analytic" && !isTRUE(use_hcv) &&
        factorization == "person_node" &&
        factor_sampling == "size_biased")) {
    cli::cli_abort(paste0(
      "Research node-box first-order geometry requires analytic, no-HCV, ",
      "person-node factorization and size-biased factor sampling."))
  }
  if (geometry == "cell_marked_node_box_first_order" &&
      !identical(as.integer(box_groups), 2L)) {
    cli::cli_abort(paste0(
      "Research cell-marked node-box geometry currently requires ",
      "{.arg box_groups = 2}."))
  }
  if (direct_cv_cell_scope == "event_horizon" &&
      !geometry %in% c("direct_person_node_cv",
                       "direct_factor_roof_person_hcv")) {
    cli::cli_abort(paste0(
      "Research {.arg direct_cv_cell_scope = \"event_horizon\"} is only ",
      "available for direct first-order or whole-person HCV geometry."))
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
  if (is.null(hcv_damping)) hcv_damping <- d_model
  if (!is.numeric(hcv_damping) || length(hcv_damping) != 1L ||
      is.na(hcv_damping) || hcv_damping <= 0) {
    cli::cli_abort("Argument {.arg hcv_damping} must be one positive number.")
  }
  if (hcv_bank_learning %in% c("exact_variance",
                               "selector_coupled_variance") &&
      !(isTRUE(use_hcv) && is.infinite(hcv_damping) &&
        identical(hcv_warmup, "hcv"))) {
    cli::cli_abort(paste0(
      "Research variance-based {.arg hcv_bank_learning} requires ",
      "active undamped HCV throughout warmup."))
  }
  if (is.infinite(hcv_damping) &&
      !(isTRUE(use_hcv) && geometry %in% c(
        "pattern_local", "norm_local_hcv_compact", "direct_person_node_hcv",
        "direct_factor_roof_person_hcv",
        "grouped_box", "node_grouped_box",
        "cell_cached_node_grouped_box",
        "lazy_node_grouped_box"))) {
    cli::cli_abort(paste0(
      "Infinite {.arg hcv_damping} (the exact undamped HCV) currently requires ",
      "{.arg use_hcv = TRUE} and {.arg geometry = \"pattern_local\", ",
      "\"norm_local_hcv_compact\", \"direct_person_node_hcv\", or ",
      "{.arg geometry = \"grouped_box\"}, \"node_grouped_box\", or ",
      "\"lazy_node_grouped_box\"."))
  }
  if (geometry %in% c("grouped_box", "node_grouped_box",
                      "cell_cached_node_grouped_box",
                      "lazy_node_grouped_box") &&
      !(backend == "analytic" && isTRUE(use_hcv) &&
        is.infinite(hcv_damping) && factorization == "person" &&
        factor_sampling == "size_biased")) {
    cli::cli_abort(paste0(
      "Grouped-box OMRF clocks currently require {.arg backend = \"analytic\"}, ",
      "{.arg use_hcv = TRUE}, {.arg hcv_damping = Inf}, ",
      "{.arg factorization = \"person\"}, and ",
      "{.arg factor_sampling = \"size_biased\"}."))
  }
  incidence <- lapply(seq_len(P), function(j) {
    incident <- which(edges[, 1L] == j | edges[, 2L] == j)
    list(
      edges = incident,
      neighbours = ifelse(edges[incident, 1L] == j,
                          edges[incident, 2L], edges[incident, 1L])
    )
  })

  # The production analytic/no-HCV routes only need one selected geometry.
  # Keep this compact path before the legacy all-geometry preparation below:
  # norm-local needs the node norm weights, while covariate-local additionally
  # needs its node-wise polynomial covariate weights.  Pattern, HCV, and box
  # routes intentionally continue through the complete preparation path.
  compact_hcv <- identical(geometry, "norm_local_hcv_compact")
  direct_hcv <- identical(geometry, "direct_person_node_hcv")
  factor_roof_person_hcv <- identical(
    geometry, "direct_factor_roof_person_hcv")
  direct_cv <- identical(geometry, "direct_person_node_cv")
  if (backend == "analytic" &&
      ((!isTRUE(use_hcv) && geometry %in% c(
        "norm_local", "direct_person_node_cv", "node_box_first_order",
        "cell_marked_node_box_first_order",
        "covariate_local",
        "covariate_local_compact")) ||
       (isTRUE(use_hcv) &&
          (compact_hcv || direct_hcv || factor_roof_person_hcv)))) {
    neighbour_norms <- NULL
    if (geometry %in% c("norm_local", "norm_local_hcv_compact",
                        "covariate_local", "covariate_local_compact")) {
      neighbour_norms <- matrix(0, nrow = P, ncol = nrow(X))
      for (node in seq_len(P)) {
        neighbours <- incidence[[node]]$neighbours
        if (length(neighbours)) {
          neighbour_norms[node, ] <- sqrt(rowSums(
            X[, neighbours, drop = FALSE]^2))
        }
      }
    }
    compact_person_node <- factorization == "person_node" &&
      geometry %in% c("norm_local", "norm_local_hcv_compact",
                      "node_box_first_order", "cell_marked_node_box_first_order")
    factor_weights <- if (compact_person_node) {
      NULL
    } else if (factorization == "person_node") {
      matrix(0, nrow = P, ncol = nrow(X) * P)
    } else {
      matrix(1, nrow = P, ncol = nrow(X))
    }
    if (factorization == "person_node" && !compact_person_node) {
      for (node in seq_len(P)) {
        columns <- (node - 1L) * nrow(X) + seq_len(nrow(X))
        factor_weights[node, columns] <- 1
      }
    }
    node_weights <- factor_weights
    norm_weights <- NULL
    if (geometry %in% c("norm_local", "norm_local_hcv_compact")) {
      powers <- if (compact_hcv) 0:3 else 0:2
      norm_weights <- matrix(0, nrow = length(powers) * P,
                             ncol = nrow(X))
      if (factorization == "person_node") {
        for (node in seq_len(P)) {
          columns <- seq_len(nrow(X))
          for (power in powers) {
            norm_weights[power * P + node, columns] <-
              neighbour_norms[node, ]^power
          }
        }
      } else {
        for (power in powers) {
          norm_weights[power * P + seq_len(P), ] <- neighbour_norms^power
        }
      }
    }
    envelope <- if (compact_person_node) {
      structure(list(type = "omrf"),
                class = c("omrf_residual_envelope", "stan_residual_envelope"))
    } else {
      stan_residual_envelope(factor_weights)
    }
    covariate_weights <- NULL
    covariate_component_blocks <- NULL
    if (geometry == "covariate_local") {
      covariate_blocks <- .omrf_covariate_blocks(X, incidence)
      covariate_weights <- do.call(rbind, covariate_blocks)
      covariate_component_blocks <- rep.int(
        seq_len(P), vapply(covariate_blocks, nrow, integer(1L)))
    }
    envelope$type <- "omrf"
    envelope$bound_type <- if (geometry %in% c(
        "node_box_first_order", "cell_marked_node_box_first_order")) {
      geometry
    } else if (direct_hcv) {
      "direct_person_node_hcv"
    } else if (factor_roof_person_hcv) {
      "direct_factor_roof_person_hcv"
    } else if (direct_cv) {
      "direct_person_node_cv"
    } else if (compact_hcv) {
      "norm_local_hcv_compact"
    } else if (geometry == "covariate_local_compact") {
      "covariate_local_compact"
    } else if (geometry == "covariate_local") {
      "covariate_local_expansion"
    } else {
      "node_local_norm_low_rank"
    }
    envelope$person_signed_roof <- person_signed_roof
    if (signed_only_person_roof)
      envelope$signed_only_person_roof <- TRUE
    if (joint_signed_group_size > 0L)
      envelope$joint_signed_group_size <- as.integer(joint_signed_group_size)
    if (shared_current_categorical)
      envelope$shared_current_categorical <- TRUE
    envelope$X <- X
    envelope$seen <- seen
    envelope$thresholds <- thresholds
    envelope$interactions <- interactions
    envelope$backend <- backend
    envelope$factorization <- factorization
    envelope$factor_sampling <- factor_sampling
    envelope$prior_backend <- prior_backend
    envelope$n_persons <- nrow(X)
    envelope$n_nodes <- P
    envelope$edge_order <- edges
    envelope$node_weights <- unname(node_weights)
    envelope$norm_weights <- unname(norm_weights)
    envelope$covariate_weights <- unname(covariate_weights)
    envelope$covariate_component_blocks <- covariate_component_blocks
    envelope$covariate_compact <- identical(
      geometry, "covariate_local_compact")
    envelope$use_hcv <- compact_hcv || direct_hcv || factor_roof_person_hcv
    envelope$hcv_after_warmup <-
      (compact_hcv || direct_hcv || factor_roof_person_hcv) &&
      identical(hcv_warmup, "first_order")
    envelope$hcv_damping <- as.numeric(hcv_damping)
    envelope$hcv_scope <- hcv_scope
    envelope$hcv_bank_learning <- hcv_bank_learning
    envelope$zigzag_sign_aware <- isTRUE(zigzag_sign_aware)
    envelope$exact_selected_refinement <-
      isTRUE(exact_selected_refinement)
    # Preserve the exact legacy immutable-envelope identity when the
    # research-only cap is disabled, so existing frozen banks remain valid.
    if (isTRUE(categorical_residual_cap)) {
      envelope$categorical_residual_cap <- TRUE
    }
    envelope$signed_bps <- isTRUE(signed_bps)
    envelope$signed_harmonic <- isTRUE(signed_harmonic)
    envelope$direct_cv_cell_scope <- direct_cv_cell_scope
    envelope$hcv_remainder_weights <- NULL
    if (geometry %in% c("node_box_first_order",
                        "cell_marked_node_box_first_order")) {
      boxes <- .omrf_node_support_boxes(X, incidence, box_groups)
      envelope$node_box_membership <- unname(boxes$membership)
      envelope$node_box_nodes <- unname(boxes$nodes)
      envelope$node_box_minima <- unname(boxes$minima)
      envelope$node_box_maxima <- unname(boxes$maxima)
      # Exact selected-factor evaluation is intrinsic to this research route.
      envelope$exact_selected_refinement <- TRUE
    }
    class(envelope) <- c("omrf_residual_envelope", "stan_residual_envelope")
    return(envelope)
  }
  row_keys <- apply(X, 1L, paste, collapse = ",")
  unique_rows <- which(!duplicated(row_keys))
  unique_weights <- numeric(length(unique_rows))
  unique_node_weights <- matrix(0, nrow = P, ncol = length(unique_rows))
  unique_hcv_remainder_weights <- matrix(0, nrow = P, ncol = length(unique_rows))
  for (key_index in seq_along(unique_rows)) {
    n <- unique_rows[key_index]
    person_design <- matrix(0, no_thresholds, d_model)
    design_row <- 1L
    for (j in seq_len(P)) {
      q <- seen[j] - 1L
      if (q == 0L) next
      B <- matrix(0, q, d_model)
      for (u in seq_len(q)) {
        B[u, threshold_starts[j] + u] <- 1
        B[u, no_thresholds + incidence[[j]]$edges] <-
          u * X[n, incidence[[j]]$neighbours]
      }
      Bnorm <- norm(B, type = "2")
      unique_node_weights[j, key_index] <- 0.5 * Bnorm^2
      # For p = softmax(z), C(p) = diag(p) - p p'. The identities
      # ||C(p)|| <= 1/2 and
      # ||C(p)-C(q)|| <= 3 ||p-q|| <= 3/2 ||z-w||
      # give the integrated categorical Taylor remainder
      # ||B'(p(Bx)-p(Ba))-B'C(p(Ba))B(x-a)||
      # <= (3/4) ||B||^3 ||x-a||^2.
      unique_hcv_remainder_weights[j, key_index] <- 0.75 * Bnorm^3
      person_design[design_row + seq_len(q) - 1L, ] <- B
      design_row <- design_row + q
    }
    unique_weights[key_index] <- 0.5 * norm(person_design, type = "2")^2
  }
  weights <- unique_weights[match(row_keys, row_keys[unique_rows])]
  node_weights <- unique_node_weights[, match(row_keys, row_keys[unique_rows]),
                                      drop = FALSE]
  hcv_remainder_weights <- unique_hcv_remainder_weights[
    , match(row_keys, row_keys[unique_rows]), drop = FALSE]
  neighbour_norms <- matrix(0, nrow = P, ncol = nrow(X))
  for (node in seq_len(P)) {
    neighbours <- incidence[[node]]$neighbours
    if (length(neighbours)) {
      neighbour_norms[node, ] <- sqrt(rowSums(
        X[, neighbours, drop = FALSE]^2))
    }
  }
  norm_weights <- matrix(0, nrow = 3L * P,
                         ncol = nrow(X))
  if (factorization == "person_node") {
    for (node in seq_len(P)) {
      columns <- seq_len(nrow(X))
      norm_weights[node, columns] <- 1
      norm_weights[P + node, columns] <- neighbour_norms[node, ]
      norm_weights[2L * P + node, columns] <- neighbour_norms[node, ]^2
    }
  } else {
    norm_weights[seq_len(P), ] <- 1
    norm_weights[P + seq_len(P), ] <- neighbour_norms
    norm_weights[2L * P + seq_len(P), ] <- neighbour_norms^2
  }
  covariate_blocks <- .omrf_covariate_blocks(X, incidence)
  # Person-node factors form a block-diagonal conceptual matrix. Store only
  # its N local columns; Julia combines these weights with the row-to-node map
  # below, avoiding P-fold structural-zero storage and alias tables.
  covariate_weights <- do.call(rbind, covariate_blocks)
  covariate_component_blocks <- rep.int(
    seq_len(P), vapply(covariate_blocks, nrow, integer(1L)))
  # The structural envelope groups factors that have the same node-local
  # neighbour configuration.  Its Julia scale callback uses the exact range
  # of the displaced categorical logits, rather than collapsing the design to
  # one spectral norm.  A factor belongs to exactly one component.
  component_nodes <- integer()
  component_covariates <- matrix(integer(), nrow = 0L, ncol = P)
  component_membership <- vector("list", P)
  for (node in seq_len(P)) {
    neighbours <- incidence[[node]]$neighbours
    keys <- if (length(neighbours)) {
      apply(X[, neighbours, drop = FALSE], 1L, paste, collapse = ",")
    } else {
      rep.int("", nrow(X))
    }
    unique_keys <- unique(keys)
    membership <- match(keys, unique_keys)
    component_membership[[node]] <- membership
    for (pattern in seq_along(unique_keys)) {
      representative <- which(membership == pattern)[1L]
      component_nodes <- c(component_nodes, node)
      component_covariates <- rbind(component_covariates, X[representative, ])
    }
  }
  component_offsets <- cumsum(c(0L, vapply(component_membership,
                                            max, integer(1L))))
  n_components <- tail(component_offsets, 1L)
  structural_groups <- if (factorization == "person_node") {
    matrix(integer(nrow(X) * P), nrow = 1L)
  } else {
    matrix(integer(P * nrow(X)), nrow = P)
  }
  for (node in seq_len(P)) {
    rows <- component_offsets[node] + component_membership[[node]]
    columns <- if (factorization == "person_node") {
      (node - 1L) * nrow(X) + seq_len(nrow(X))
    } else {
      seq_len(nrow(X))
    }
    if (factorization == "person_node") {
      structural_groups[1L, columns] <- rows
    } else {
      structural_groups[node, columns] <- rows
    }
  }
  if (factorization == "person_node") {
    factor_weights <- matrix(0, nrow = P, ncol = nrow(X) * P)
    factor_hcv_weights <- matrix(0, nrow = P, ncol = nrow(X) * P)
    for (node in seq_len(P)) {
      columns <- (node - 1L) * nrow(X) + seq_len(nrow(X))
      factor_weights[node, columns] <- node_weights[node, ]
      factor_hcv_weights[node, columns] <- hcv_remainder_weights[node, ]
    }
    envelope <- stan_residual_envelope(factor_weights)
    node_weights <- factor_weights
    hcv_remainder_weights <- factor_hcv_weights
  } else {
    envelope <- stan_residual_envelope(weights)
  }

  # Deterministic farthest-first seeds avoid changing the caller's RNG state.
  # The subsequent Lloyd iterations are deterministic because explicit
  # distinct centers are supplied. Boxes need not approximate a probability
  # model: any fixed partition is exact because its coordinate extrema
  # certify every person assigned to it.
  box_membership <- rep.int(1L, nrow(X))
  box_minima <- matrix(apply(X, 2L, min), nrow = 1L)
  box_maxima <- matrix(apply(X, 2L, max), nrow = 1L)
  partition_boxes <- function(features, group_limit) {
    if (!ncol(features)) {
      return(list(membership = rep.int(1L, nrow(features)), groups = 1L))
    }
    keys <- apply(features, 1L, paste, collapse = ",")
    distinct <- features[!duplicated(keys), , drop = FALSE]
    group_count <- min(group_limit, nrow(distinct))
    center_rows <- integer(group_count)
    center_rows[[1L]] <- which.max(rowSums(
      sweep(distinct, 2L, colMeans(distinct), "-")^2))
    minimum_distance <- rowSums(
      sweep(distinct, 2L, distinct[center_rows[[1L]], ], "-")^2)
    if (group_count > 1L) {
      for (group in 2:group_count) {
        center_rows[[group]] <- which.max(minimum_distance)
        distance <- rowSums(
          sweep(distinct, 2L, distinct[center_rows[[group]], ], "-")^2)
        minimum_distance <- pmin(minimum_distance, distance)
      }
    }
    membership <- if (group_count == 1L) rep.int(1L, nrow(features)) else {
      stats::kmeans(
        features, centers = distinct[center_rows, , drop = FALSE],
        iter.max = 100L)$cluster
    }
    list(membership = as.integer(membership), groups = group_count)
  }
  partition_support_tree <- function(features, group_limit, root_limit) {
    if (!ncol(features)) {
      return(list(
        membership = rep.int(1L, nrow(features)), groups = 1L,
        boxes = list(seq_len(nrow(features))), parents = 0L,
        root_nodes = 1L, leaves = 1L))
    }
    boxes <- list(seq_len(nrow(features)))
    parents <- 0L
    active <- 1L
    root_nodes <- if (root_limit == 1L) 1L else integer()
    while (length(active) < group_limit) {
      candidates <- lapply(active, function(index) {
        rows <- boxes[[index]]
        ranges <- apply(features[rows, , drop = FALSE], 2L, function(x) {
          diff(range(x))
        })
        list(index = index, score = sum(ranges) * length(rows), ranges = ranges)
      })
      splittable <- which(vapply(
        candidates, function(candidate) any(candidate$ranges > 0), logical(1L)))
      if (!length(splittable)) break
      selected_candidate <- splittable[[which.max(vapply(
        candidates[splittable], `[[`, numeric(1L), "score"))]]
      selected <- candidates[[selected_candidate]]$index
      rows <- boxes[[selected]]
      dimension <- which.max(candidates[[selected_candidate]]$ranges)
      ordered <- rows[order(features[rows, dimension], rows)]
      changes <- which(diff(features[ordered, dimension]) > 0)
      split <- changes[[which.min(abs(changes - length(rows) / 2))]]
      left <- ordered[seq_len(split)]
      right <- ordered[seq.int(split + 1L, length(ordered))]
      if (!length(left) || !length(right)) {
        stop("Internal node support-box split failed: selected=", selected,
             ", rows=", length(rows), ", dimension=", dimension,
             ", changes=", paste(changes, collapse = ","),
             ", split=", split, ", boxes=", length(boxes), call. = FALSE)
      }
      boxes <- c(boxes, list(left, right))
      parents <- c(parents, selected, selected)
      active <- c(active[active != selected], length(boxes) - 1L,
                  length(boxes))
      if (!length(root_nodes) && length(active) >= root_limit) {
        root_nodes <- active
      }
    }
    if (!length(root_nodes)) root_nodes <- active
    membership <- integer(nrow(features))
    for (group in active) membership[boxes[[group]]] <- group
    list(membership = membership, groups = length(boxes), boxes = boxes,
         parents = parents, root_nodes = root_nodes, leaves = active)
  }
  partition_support_boxes <- function(features, group_limit) {
    tree <- partition_support_tree(features, group_limit, group_limit)
    leaf_map <- integer(tree$groups)
    leaf_map[tree$leaves] <- seq_along(tree$leaves)
    list(membership = leaf_map[tree$membership],
         groups = length(tree$leaves))
  }
  node_box_membership <- matrix(1L, nrow = P, ncol = nrow(X))
  node_box_nodes <- seq_len(P)
  node_box_minima <- matrix(0, nrow = P, ncol = P)
  node_box_maxima <- matrix(0, nrow = P, ncol = P)
  node_box_parents <- integer()
  node_box_roots <- integer()
  node_box_leaves <- integer()
  if (identical(geometry, "grouped_box")) {
    partition <- partition_boxes(X, box_groups)
    group_count <- partition$groups
    box_membership <- partition$membership
    box_minima <- do.call(rbind, lapply(seq_len(group_count), function(group) {
      apply(X[box_membership == group, , drop = FALSE], 2L, min)
    }))
    box_maxima <- do.call(rbind, lapply(seq_len(group_count), function(group) {
      apply(X[box_membership == group, , drop = FALSE], 2L, max)
    }))
  } else if (geometry %in% c("node_grouped_box",
                             "cell_cached_node_grouped_box",
                             "lazy_node_grouped_box")) {
    memberships <- vector("list", P)
    minima <- list()
    maxima <- list()
    node_box_nodes <- integer()
    offset <- 0L
    for (node in seq_len(P)) {
      neighbours <- incidence[[node]]$neighbours
      features <- X[, neighbours, drop = FALSE]
      lazy <- identical(geometry, "lazy_node_grouped_box")
      partition <- if (lazy) {
        partition_support_tree(features, box_groups, box_root_groups)
      } else {
        partition_support_boxes(features, box_groups)
      }
      local_membership <- partition$membership
      memberships[[node]] <- offset + local_membership
      for (group in seq_len(partition$groups)) {
        rows <- if (lazy) seq_len(nrow(X)) %in% partition$boxes[[group]] else
          local_membership == group
        lower <- upper <- numeric(P)
        if (length(neighbours)) {
          lower[neighbours] <- apply(
            features[rows, , drop = FALSE], 2L, min)
          upper[neighbours] <- apply(
            features[rows, , drop = FALSE], 2L, max)
        }
        minima[[offset + group]] <- lower
        maxima[[offset + group]] <- upper
      }
      node_box_nodes <- c(
        node_box_nodes, rep.int(node, partition$groups))
      if (lazy) {
        node_box_parents <- c(node_box_parents,
          ifelse(partition$parents == 0L, 0L, offset + partition$parents))
        node_box_roots <- c(node_box_roots, offset + partition$root_nodes)
        node_box_leaves <- c(node_box_leaves, offset + partition$leaves)
      }
      offset <- offset + partition$groups
    }
    node_box_membership <- do.call(rbind, memberships)
    node_box_minima <- do.call(rbind, minima)
    node_box_maxima <- do.call(rbind, maxima)
  }
  envelope$type <- "omrf"
  envelope$bound_type <- if (backend == "analytic" && use_hcv &&
                             geometry %in% c("node_grouped_box",
                                             "cell_cached_node_grouped_box",
                                             "lazy_node_grouped_box") &&
                             is.infinite(hcv_damping)) {
    if (geometry == "lazy_node_grouped_box") {
      "lazy_node_grouped_person_hcv_box"
    } else if (geometry == "cell_cached_node_grouped_box") {
      "cell_cached_node_grouped_person_hcv_box"
    } else {
      "node_grouped_person_hcv_box"
    }
  } else if (backend == "analytic" && use_hcv &&
                             geometry == "grouped_box" &&
                             is.infinite(hcv_damping)) {
    "grouped_person_hcv_box"
  } else if (backend == "analytic" && use_hcv &&
                             geometry == "pattern_local" &&
                             is.infinite(hcv_damping)) {
    "pattern_local_hcv_range"
  } else if (backend == "analytic" && !use_hcv &&
                             geometry == "pattern_local") {
    "pattern_local_range"
  } else if (backend == "analytic" && !use_hcv &&
             geometry == "covariate_local") {
    "covariate_local_expansion"
  } else if (backend == "analytic" && !use_hcv) {
    "node_local_norm_low_rank"
  } else {
    "stacked_person_spectral"
  }
  envelope$X <- X
  envelope$seen <- seen
  envelope$thresholds <- thresholds
  envelope$interactions <- interactions
  envelope$backend <- backend
  envelope$factorization <- factorization
  envelope$factor_sampling <- factor_sampling
  envelope$prior_backend <- prior_backend
  envelope$n_persons <- nrow(X)
  envelope$n_nodes <- P
  envelope$edge_order <- edges
  envelope$node_weights <- unname(node_weights)
  envelope$norm_weights <- unname(norm_weights)
  envelope$covariate_weights <- unname(covariate_weights)
  envelope$covariate_component_blocks <- covariate_component_blocks
  envelope$structural_groups <- unname(structural_groups)
  envelope$n_structural_components <- n_components
  envelope$component_nodes <- component_nodes
  envelope$component_covariates <- unname(component_covariates)
  envelope$use_hcv <- isTRUE(use_hcv)
  envelope$hcv_after_warmup <- isTRUE(use_hcv) &&
    identical(hcv_warmup, "first_order")
  envelope$hcv_damping <- as.numeric(hcv_damping)
  envelope$hcv_scope <- hcv_scope
  envelope$hcv_bank_learning <- hcv_bank_learning
  envelope$zigzag_sign_aware <- isTRUE(zigzag_sign_aware)
  envelope$exact_selected_refinement <-
    isTRUE(exact_selected_refinement)
  envelope$signed_bps <- isTRUE(signed_bps)
  envelope$signed_harmonic <- isTRUE(signed_harmonic)
  envelope$hcv_remainder_weights <- unname(hcv_remainder_weights)
  envelope$box_membership <- as.integer(box_membership)
  envelope$box_minima <- unname(box_minima)
  envelope$box_maxima <- unname(box_maxima)
  envelope$node_box_membership <- unname(node_box_membership)
  envelope$node_box_nodes <- as.integer(node_box_nodes)
  envelope$node_box_minima <- unname(node_box_minima)
  envelope$node_box_maxima <- unname(node_box_maxima)
  envelope$node_box_parents <- as.integer(node_box_parents)
  envelope$node_box_roots <- as.integer(node_box_roots)
  envelope$node_box_leaves <- as.integer(node_box_leaves)
  class(envelope) <- c("omrf_residual_envelope", "stan_residual_envelope")
  envelope
}

#' Native full-gradient backend for an analytic-prior OMRF
#'
#' Constructs a native Julia backend for the complete, non-factorized OMRF
#' negative gradient and Hessian-vector product. It is intended for full
#' Adaptive Boomerang, Boomerang, BPS, and Zig-Zag samplers; it does not change
#' the sampler into factorized dynamics. The Stan model is still used for
#' unconstrained parameter metadata and output conversion.
#'
#' The model data passed to [pdmp_sample_from_stanmodel()] must contain the
#' threshold-prior parameters and either the independent interaction scale or
#' the node-scale Gaussian parameters. For data objects carrying both prior
#' families, pass the prior-only data through `prior_standata` so the analytic
#' node-scale branch is selected explicitly. Threshold, interaction, and
#' node-scale blocks must cover the complete unconstrained state.
#'
#' @inheritParams omrf_residual_envelope
#' @param prior_standata Optional prior-only OMRF data. Supplying this is
#'   required when the full data object carries both independent and dependent
#'   prior fields and the node-scale Gaussian branch is intended.
#' @return An OMRF full-gradient backend specification for the
#'   `full_gradient` argument of [pdmp_sample_from_stanmodel()].
#' @export
omrf_full_gradient <- function(X, seen, thresholds, interactions,
                               prior_standata = NULL) {
  spec <- omrf_residual_envelope(
    X, seen, thresholds, interactions,
    backend = "analytic", use_hcv = FALSE,
    factorization = "person", factor_sampling = "size_biased",
    prior_backend = "analytic")
  if (!is.null(prior_standata)) {
    spec <- .attach_omrf_analytic_prior(
      spec, prior_standata, "prior_standata")
  }
  structure(list(spec = spec), class = "omrf_full_gradient")
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
#' @rdname stan_subsampling
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

#' Diagnose a subsampling custom-Stan gradient and envelope
#'
#' Evaluates the persistent full, prior-only, and selected Stan modes at one
#' position without starting a sampler. The selected likelihood gradient is
#' unscaled; `selected_scale` reports the `N / m` factor applied by Julia when
#' forming the subsampling gradient. The residual event-rate diagnostic uses the
#' same typed envelope and trajectory geometry as sampling.
#'
#' @param path_to_stanmodel Stan source or compiled-library path containing the
#'   package subset hook.
#' @param standata Full-likelihood Stan data list or JSON path.
#' @param subsampling A specification from [stan_subsampling()].
#' @param position Finite unconstrained parameter vector.
#' @param subset One-based observation indices, with exactly the configured
#'   subsample size and no duplicates.
#' @param velocity Finite diagnostic velocity. Defaults to a vector of ones.
#' @param flow Supported subsampling flow used for the rate diagnostic.
#' @param flow_mean Optional flow mean.
#' @param flow_cov Optional flow covariance.
#' @param phase Diagnostic sampling phase. For staged HCV specifications,
#'   `"main"` activates HCV and `"warmup"` uses the first-order fallback.
#' @return A list of gradient-closure, residual-rate, envelope, and construction
#'   diagnostics.
#' @export
stan_subsampling_diagnostics <- function(path_to_stanmodel, standata,
                                    subsampling, position, subset,
                                    velocity = NULL,
                                    flow = c("ZigZag", "BouncyParticle",
                                             "AdaptiveBoomerang"),
                                    flow_mean = NULL, flow_cov = NULL,
                                    phase = c("main", "warmup")) {
  validate_type(path_to_stanmodel, type = "character", n = 1)
  if (!inherits(subsampling, "stan_subsampling")) {
    cli::cli_abort("Argument {.arg subsampling} must be created by {.fn stan_subsampling}.")
  }
  if (!is.numeric(position) || !length(position) || any(!is.finite(position))) {
    cli::cli_abort("Argument {.arg position} must be a nonempty finite numeric vector.")
  }
  if (!rlang::is_integerish(subset, n = subsampling$subsample_size,
                            finite = TRUE) ||
      any(subset < 1L | subset > subsampling$n_observations) ||
      anyDuplicated(subset)) {
    cli::cli_abort("Argument {.arg subset} must contain exactly m distinct indices in 1:N.")
  }
  if (is.null(velocity)) velocity <- rep(1, length(position))
  if (!is.numeric(velocity) || length(velocity) != length(position) ||
      any(!is.finite(velocity))) {
    cli::cli_abort("Argument {.arg velocity} must be finite and match {.arg position}.")
  }
  flow <- match.arg(flow)
  phase <- match.arg(phase)
  if (is.null(flow_mean)) flow_mean <- rep(0, length(position))
  if (is.null(flow_cov)) flow_cov <- diag(length(position))

  write_data <- function(data) {
    if (!is.list(data)) {
      validate_type(data, type = "character", n = 1)
      return(normalizePath(data, mustWork = TRUE))
    }
    path <- tempfile(fileext = ".json")
    write_stan_json(data, path)
    path
  }
  full_data <- write_data(standata)
  prior_data <- write_data(subsampling$prior_standata)
  temporary_data <- c(
    if (is.list(standata)) full_data else character(),
    if (is.list(subsampling$prior_standata)) prior_data else character()
  )
  if (length(temporary_data)) on.exit(unlink(temporary_data), add = TRUE)

  check_for_julia_setup()
  model_path <- normalizePath(path_to_stanmodel, mustWork = TRUE)
  if (grepl("\\.stan$", model_path)) {
    model_path <- compile_pdmp_stan_model(model_path)
  }
  JuliaCall::julia_assign("_diagnostic_model", model_path)
  JuliaCall::julia_assign("_diagnostic_full_data", full_data)
  JuliaCall::julia_assign("_diagnostic_prior_data", prior_data)
  JuliaCall::julia_assign("_diagnostic_subsampling", subsampling)
  JuliaCall::julia_assign("_diagnostic_position", as.numeric(position))
  JuliaCall::julia_assign("_diagnostic_subset", as.integer(subset))
  JuliaCall::julia_assign("_diagnostic_velocity", as.numeric(velocity))
  JuliaCall::julia_assign("_diagnostic_flow", flow)
  JuliaCall::julia_assign("_diagnostic_flow_mean", as.numeric(flow_mean))
  JuliaCall::julia_assign("_diagnostic_flow_cov", flow_cov)
  JuliaCall::julia_assign("_diagnostic_phase", phase)
  result <- .pdmpsamplers_julia_eval(
    "PDMPSamplersRBridge.r_stan_subsampling_diagnostics(
      _diagnostic_model, _diagnostic_full_data, _diagnostic_prior_data,
      _diagnostic_subsampling, _diagnostic_position, _diagnostic_subset,
      _diagnostic_velocity, _diagnostic_flow, _diagnostic_flow_mean,
      _diagnostic_flow_cov, _diagnostic_phase
    );"
  )
  if (is.environment(result)) result <- as.list(result)
  result
}
