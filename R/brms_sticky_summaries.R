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
    cli::cli_warn(
        c(
            "The median probability model can behave poorly under moderate or strong collinearity.",
            "i" = "Coefficient-wise inclusion and term-wise predictive model reduction are not the same task."
        ),
        class = "mpm_collinearity_warning"
    )
    ip <- inclusion_prob(x, ...)
    names(ip)[ip > threshold]
}

# ──────────────────────────────────────────────────────────────────────────────
# model_averaged_mean
# ──────────────────────────────────────────────────────────────────────────────

#' Model-averaged posterior means for stickable coefficients
#'
#' Returns \eqn{E[\beta_k \mid y]} for each stickable population-level
#' coefficient. Because the sticky PDMP trajectory includes time at zero
#' (the spike), this is the spike-and-slab model-averaged mean, not the
#' conditional mean given inclusion.
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
        ggplot2::geom_hline(yintercept = threshold, linetype = "dashed", colour = "grey40") +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = NULL, y = "Inclusion probability")
}
