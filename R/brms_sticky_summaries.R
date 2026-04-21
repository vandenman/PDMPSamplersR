# Post-processing summaries for sticky brms fits.
#
# Implements coefficient-wise inclusion probability accessors,
# the median probability model, model-averaged summaries, an
# inclusion plot for brmsfit objects fitted with
# brm_pdmp(..., sticky = TRUE).

.check_brmsfit_sticky <- function(x) {
    if (!inherits(x, "brmsfit"))
        cli::cli_abort("{.arg x} must be a {.cls brmsfit} object.")
    sticky_meta <- attr(x, "sticky")
    if (is.null(sticky_meta))
        cli::cli_abort(c(
            "{.arg x} does not have sticky variable selection metadata.",
            "i" = "Run {.fn brm_pdmp} with {.code sticky = TRUE} to enable variable selection."
        ))
    sticky_meta
}

# ──────────────────────────────────────────────────────────────────────────────
# inclusion_prob
# ──────────────────────────────────────────────────────────────────────────────

#' Coefficient-wise posterior inclusion probabilities
#'
#' Extracts per-coefficient inclusion probabilities from a sticky
#' `brmsfit` fitted with [brm_pdmp()] and `sticky = TRUE`. Each
#' probability estimates \eqn{P(\gamma_k = 1 \mid y)}, the posterior
#' probability that coefficient \eqn{\beta_k} is non-zero under the
#' user-specified spike-and-slab prior.
#'
#' @param x A `brmsfit` fitted with `sticky = TRUE`.
#' @param chain Integer, which chain to extract, or `NULL` to average
#'   across all chains (default: `NULL`).
#' @param ... Ignored.
#'
#' @return A named numeric vector in \eqn{[0, 1]}, one entry per
#'   stickable population-level coefficient.
#' @export
inclusion_prob <- function(x, ...) UseMethod("inclusion_prob")

#' @rdname inclusion_prob
#' @export
inclusion_prob.default <- function(x, ...) {
    cli::cli_abort("{.arg x} must be a {.cls brmsfit} object with sticky metadata.")
}

#' @rdname inclusion_prob
#' @export
inclusion_prob.brmsfit <- function(x, chain = NULL, ...) {
    sticky_meta <- .check_brmsfit_sticky(x)
    incl_list <- sticky_meta$inclusion_probs
    if (is.null(chain)) {
        m <- do.call(rbind, incl_list)
        colMeans(m)
    } else {
        chain_nm <- paste0("chain", chain)
        if (!chain_nm %in% names(incl_list))
            cli::cli_abort("Chain {.val {chain}} not found in the sticky result.")
        incl_list[[chain_nm]]
    }
}

# ──────────────────────────────────────────────────────────────────────────────
# inclusion_bf
# ──────────────────────────────────────────────────────────────────────────────

#' Inclusion Bayes factors for stickable coefficients
#'
#' Computes the inclusion Bayes factor for each stickable coefficient:
#' \deqn{BF_{\text{incl},k} = \frac{P(\gamma_k = 1 \mid y)}{P(\gamma_k = 0 \mid y)}
#'       \cdot \frac{P(\gamma_k = 0)}{P(\gamma_k = 1)}}
#'
#' This answers "how much did the data update the prior odds in favour of
#' inclusion?" Values greater than 1 indicate evidence for inclusion; values
#' less than 1 indicate evidence against. Unlike inclusion probabilities alone,
#' the Bayes factor separates the likelihood signal from the prior.
#'
#' @param x A `brmsfit` fitted with `sticky = TRUE` using a version of
#'   [brm_pdmp()] that stores the model prior in the sticky metadata.
#' @param ... Passed to [inclusion_prob()].
#'
#' @return A named numeric vector of inclusion Bayes factors, one per
#'   stickable coefficient.
#' @export
inclusion_bf <- function(x, ...) UseMethod("inclusion_bf")

#' @rdname inclusion_bf
#' @export
inclusion_bf.default <- function(x, ...) {
    cli::cli_abort("{.arg x} must be a {.cls brmsfit} object with sticky metadata.")
}

#' @rdname inclusion_bf
#' @export
inclusion_bf.brmsfit <- function(x, ...) {
    sticky_meta <- .check_brmsfit_sticky(x)
    model_prior <- sticky_meta$model_prior
    if (is.null(model_prior))
        cli::cli_abort(c(
            "No model prior found in sticky metadata.",
            "i" = "Refit with the current version of {.fn brm_pdmp} to include model prior information."
        ))
    ip <- inclusion_prob(x, ...)
    stickable_idx <- which(sticky_meta$can_stick)
    if (is.bernoulli(model_prior)) {
        pi_prior <- model_prior$prob[stickable_idx]
    } else {
        # betabernoulli: marginal inclusion probability = a / (a + b)
        pi_prior <- rep(model_prior$a / (model_prior$a + model_prior$b), length(ip))
    }
    names(pi_prior) <- names(ip)
    (ip / (1 - ip)) / (pi_prior / (1 - pi_prior))
}

# ──────────────────────────────────────────────────────────────────────────────
# median_probability_model
# ──────────────────────────────────────────────────────────────────────────────

#' Median probability model
#'
#' Returns the names of population-level coefficients with posterior
#' inclusion probability above `threshold` (default 0.5). This is the
#' *median probability model*: the model that includes each variable for
#' which the posterior probability of inclusion exceeds one half.
#'
#' @param x A `brmsfit` fitted with `sticky = TRUE`.
#' @param threshold Numeric scalar in \eqn{[0, 1]} (default: 0.5).
#' @param ... Passed to [inclusion_prob()].
#'
#' @return A character vector of coefficient names.
#'
#' @section Collinearity caveat:
#' The median probability model can behave poorly under moderate or
#' strong collinearity. Coefficient-wise inclusion and term-wise
#' predictive reduction are not the same task; see the package
#' documentation for guidance.
#'
#' @export
median_probability_model <- function(x, ...) UseMethod("median_probability_model")

#' @rdname median_probability_model
#' @export
median_probability_model.brmsfit <- function(x, threshold = 0.5, ...) {
    if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0 || threshold > 1)
        cli::cli_abort("{.arg threshold} must be a numeric scalar in [0, 1].")
    ip <- inclusion_prob(x, ...)
    .check_mpm_collinearity(x, ip)
    names(ip)[ip > threshold]
}

.check_mpm_collinearity <- function(x, ip, cor_threshold = 0.3) {
    stickable_nms <- paste0("b_", sub("^b\\.", "", names(ip)))
    draws_mat <- tryCatch(
        as.matrix(brms::as_draws_df(x)),
        error = function(e) NULL
    )
    if (is.null(draws_mat)) return(invisible(NULL))
    available_nms <- stickable_nms[stickable_nms %in% colnames(draws_mat)]
    if (length(available_nms) < 2L) return(invisible(NULL))
    cor_mat <- abs(stats::cor(draws_mat[, available_nms, drop = FALSE]))
    diag(cor_mat) <- 0
    if (max(cor_mat) > cor_threshold)
        cli::cli_warn(
            c(
                "The median probability model can behave poorly under moderate or strong collinearity.",
                "i" = "Coefficient-wise inclusion and term-wise predictive model reduction are not the same task."
            ),
            class = "mpm_collinearity_warning"
        )
    invisible(NULL)
}

# ──────────────────────────────────────────────────────────────────────────────
# model_averaged_mean
# ──────────────────────────────────────────────────────────────────────────────

#' Model-averaged posterior means for stickable coefficients
#'
#' Returns \eqn{E[\beta_k \mid y]} for each stickable population-level
#' coefficient by averaging the posterior draws produced by the sticky PDMP
#' trajectory. Because the trajectory spends time at zero (the spike) in
#' proportion to the posterior spike mass, this equals the spike-and-slab
#' model-averaged mean provided the sampler has converged and the zero-dwell
#' time is correctly attributed. Use [inclusion_prob()] to inspect the
#' corresponding inclusion probabilities.
#'
#' @param x A `brmsfit` fitted with `sticky = TRUE`.
#' @param ... Ignored.
#'
#' @return A named numeric vector of posterior means.
#' @importFrom rlang .data
#' @export
model_averaged_mean <- function(x, ...) UseMethod("model_averaged_mean")

#' @rdname model_averaged_mean
#' @export
model_averaged_mean.brmsfit <- function(x, ...) {
    .check_brmsfit_sticky(x)
    if (!requireNamespace("brms", quietly = TRUE))
        cli::cli_abort("Package {.pkg brms} is required for {.fn model_averaged_mean}.")
    ip_names <- names(inclusion_prob(x))
    coef_names <- sub("^b\\.", "", ip_names)
    fe <- brms::fixef(x)
    fe[rownames(fe) %in% coef_names, "Estimate"]
}

# ──────────────────────────────────────────────────────────────────────────────
# check_inclusion_stability
# ──────────────────────────────────────────────────────────────────────────────

#' Per-chain stability check for inclusion probability estimates
#'
#' Compares inclusion probability estimates across chains to help assess
#' convergence of the sticky sampler. Returns a data frame with cross-chain
#' summary statistics for each stickable coefficient.
#'
#' @param x A `brmsfit` fitted with `sticky = TRUE`.
#' @param ... Ignored.
#'
#' @return A data frame with columns `coefficient`, `mean`, `sd`, `min`,
#'   `max`, and `range` summarising cross-chain variability.
#' @export
check_inclusion_stability <- function(x, ...) UseMethod("check_inclusion_stability")

#' @rdname check_inclusion_stability
#' @export
check_inclusion_stability.default <- function(x, ...) {
    cli::cli_abort("{.arg x} must be a {.cls brmsfit} object with sticky metadata.")
}

#' @rdname check_inclusion_stability
#' @export
check_inclusion_stability.brmsfit <- function(x, ...) {
    sticky_meta <- .check_brmsfit_sticky(x)
    incl_list <- sticky_meta$inclusion_probs
    if (length(incl_list) < 2L)
        cli::cli_inform("Only one chain available; cross-chain stability cannot be assessed.")
    m <- do.call(rbind, incl_list)
    rownames(m) <- paste0("chain", seq_len(nrow(m)))
    data.frame(
        coefficient = colnames(m),
        mean        = colMeans(m),
        sd          = apply(m, 2, stats::sd),
        min         = apply(m, 2, min),
        max         = apply(m, 2, max),
        range       = apply(m, 2, function(v) diff(range(v))),
        row.names   = NULL
    )
}

# ──────────────────────────────────────────────────────────────────────────────
# plot_inclusion
# ──────────────────────────────────────────────────────────────────────────────

#' Plot posterior inclusion probabilities
#'
#' Produces a horizontal bar chart of coefficient-wise inclusion
#' probabilities from a sticky `brmsfit`, with a vertical reference
#' line at `threshold`.
#'
#' @param x A `brmsfit` fitted with `sticky = TRUE`.
#' @param threshold Numeric scalar in \eqn{[0, 1]}, the reference line
#'   position (default: 0.5).
#' @param ... Passed to [inclusion_prob()].
#'
#' @return A `ggplot` object.
#' @export
plot_inclusion <- function(x, ...) UseMethod("plot_inclusion")

#' @rdname plot_inclusion
#' @export
plot_inclusion.brmsfit <- function(x, threshold = 0.5, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
        cli::cli_abort("Package {.pkg ggplot2} is required for {.fn plot_inclusion}.")
    ip <- inclusion_prob(x, ...)
    df <- data.frame(
        coefficient = factor(names(ip), levels = names(ip)[order(ip)]),
        inclusion_prob = as.numeric(ip)
    )
    ggplot2::ggplot(df, ggplot2::aes(x = rlang::.data$coefficient, y = rlang::.data$inclusion_prob)) +
        ggplot2::geom_col() +
        ggplot2::geom_text(
            ggplot2::aes(label = sprintf("%.2f", rlang::.data$inclusion_prob)),
            hjust = -0.1, size = 3
        ) +
        ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", colour = "grey40") +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(limits = c(0, 1.15)) +
        ggplot2::labs(x = NULL, y = "Inclusion probability")
}

# ──────────────────────────────────────────────────────────────────────────────
# print.sticky_brmsfit
# ──────────────────────────────────────────────────────────────────────────────

#' Print method for sticky brmsfit objects
#'
#' Delegates to the standard brms print method and appends a one-line note
#' about sticky variable selection results.
#'
#' @param x A `sticky_brmsfit` object returned by [brm_pdmp()] with
#'   `sticky = TRUE`.
#' @param digits Number of decimal digits to print (default: 2).
#' @param short Logical; whether to use the short summary format
#'   (default: value of `getOption("brms.short_summary", FALSE)`).
#' @param ... Passed to the brms print method.
#'
#' @return `x`, invisibly.
#' @export
print.sticky_brmsfit <- function(x, digits = 2,
                                 short = getOption("brms.short_summary", FALSE), ...) {
    NextMethod()
    cli::cli_inform("Sticky variable selection was performed; call {.fn inclusion_prob} to inspect results.")
    invisible(x)
}
