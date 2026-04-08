# Sticky variable selection helpers for brm_pdmp()
#
# Maps brms formula terms to unconstrained BridgeStan coordinates and
# derives per-coefficient kappa values from brms prior specifications.

# ──────────────────────────────────────────────────────────────────────────────
# can_stick mapping
# ──────────────────────────────────────────────────────────────────────────────

#' Build logical can_stick vector from brms metadata
#'
#' Determines which unconstrained BridgeStan coordinates correspond to
#' supported population-level coefficients (excluding intercept, scale,
#' and shape parameters).
#'
#' @param unc_names Character vector from `BridgeStan::param_unc_names()`.
#' @param user_can_stick Optional logical vector (same length as supported
#'   non-intercept `b` coefficients) to override the default all-TRUE.
#'
#' @return A logical vector of length `length(unc_names)`.
#'
#' @keywords internal
map_can_stick <- function(unc_names, user_can_stick = NULL) {
    b_indices <- grep("^b\\.", unc_names)
    intercept_indices <- grep("^b\\..*[Ii]ntercept", unc_names)
    stickable_indices <- setdiff(b_indices, intercept_indices)

    if (!is.null(user_can_stick)) {
        if (!is.logical(user_can_stick))
            cli::cli_abort("{.arg can_stick} must be a logical vector.")
        if (length(user_can_stick) != length(stickable_indices))
            cli::cli_abort(c(
                "{.arg can_stick} must have length {length(stickable_indices)} (the number of supported non-intercept population-level coefficients).",
                "i" = "Stickable coefficients: {.val {unc_names[stickable_indices]}}.",
                "x" = "Got length {length(user_can_stick)}."
            ))
        out <- rep(FALSE, length(unc_names))
        out[stickable_indices] <- user_can_stick
    } else {
        out <- rep(FALSE, length(unc_names))
        out[stickable_indices] <- TRUE
    }
    out
}

#' Return the names of stickable coefficients
#'
#' @param unc_names Character vector from `BridgeStan::param_unc_names()`.
#' @return Character vector of stickable unconstrained parameter names.
#' @keywords internal
stickable_coef_names <- function(unc_names) {
    b_indices <- grep("^b\\.", unc_names)
    intercept_indices <- grep("^b\\..*[Ii]ntercept", unc_names)
    stickable_indices <- setdiff(b_indices, intercept_indices)
    unc_names[stickable_indices]
}

# ──────────────────────────────────────────────────────────────────────────────
# Automatic parameter_prior derivation
# ──────────────────────────────────────────────────────────────────────────────

#' Derive parameter_prior from brms prior specification
#'
#' Evaluates the slab density at zero for each unconstrained coordinate.
#' Only `normal(0, s)` and `student_t(df, 0, s)` priors on `b` class
#' coefficients are supported. All other cases error out.
#'
#' @param prior A `brmsprior` data frame (from `brms::prior_summary()` or
#'   the prior slot of a brmsfit).
#' @param unc_names Character vector from `BridgeStan::param_unc_names()`.
#' @param can_stick Logical vector of length `d` (from `map_can_stick()`).
#'
#' @return Numeric vector of length `d` with slab densities at zero for
#'   stickable coordinates and 1.0 elsewhere.
#'
#' @keywords internal
derive_parameter_prior <- function(prior, unc_names, can_stick) {
    d <- length(unc_names)
    param_prior <- rep(1.0, d)
    stickable_idx <- which(can_stick)
    if (length(stickable_idx) == 0L) return(param_prior)

    stickable_names <- unc_names[stickable_idx]

    for (i in seq_along(stickable_idx)) {
        coef_name <- stickable_names[i]
        coef_short <- sub("^b\\.", "", coef_name)
        density_at_zero <- .prior_density_at_zero(prior, coef_short)
        param_prior[stickable_idx[i]] <- density_at_zero
    }

    param_prior
}

#' Evaluate the slab density at zero for a single b coefficient
#'
#' Looks up the prior for `coef` first as a coefficient-specific prior,
#' then falls back to the class-level `b` prior. Supports `normal(0, s)`
#' and `student_t(df, 0, s)` only.
#'
#' @param prior A `brmsprior` data frame.
#' @param coef Character, the brms coefficient name (without `b.` prefix).
#'
#' @return Numeric scalar, the slab density evaluated at zero.
#' @keywords internal
.prior_density_at_zero <- function(prior, coef) {
    # First try coefficient-specific prior
    row <- prior[prior$class == "b" & prior$coef == coef & nchar(prior$prior) > 0, ]
    if (nrow(row) == 0L) {
        # Fall back to class-level prior
        row <- prior[prior$class == "b" & prior$coef == "" & nchar(prior$prior) > 0, ]
    }
    if (nrow(row) == 0L) {
        cli::cli_abort(c(
            "Cannot automatically derive {.arg parameter_prior} for coefficient {.val {coef}}.",
            "i" = "No explicit prior found on class {.val b}.",
            "i" = "Supply {.arg parameter_prior} explicitly."
        ))
    }

    prior_str <- row$prior[1L]
    .parse_and_eval_prior_at_zero(prior_str, coef)
}

#' Parse a brms prior string and evaluate the density at zero
#'
#' @param prior_str Character, e.g. "normal(0, 2.5)" or "student_t(3, 0, 2.5)".
#' @param coef Character, coefficient name for error messages.
#' @return Numeric scalar.
#' @keywords internal
.parse_and_eval_prior_at_zero <- function(prior_str, coef) {
    normal_match <- regmatches(
        prior_str,
        regexec("^normal\\(\\s*([^,]+)\\s*,\\s*([^)]+)\\s*\\)$", prior_str)
    )[[1L]]
    if (length(normal_match) == 3L) {
        mu <- as.numeric(normal_match[2L])
        sigma <- as.numeric(normal_match[3L])
        if (!isTRUE(mu == 0))
            cli::cli_abort(c(
                "Automatic {.arg parameter_prior} derivation requires a zero-centered slab prior for coefficient {.val {coef}}.",
                "i" = "Found {.val {prior_str}} (mean = {mu}).",
                "i" = "Supply {.arg parameter_prior} explicitly."
            ))
        return(stats::dnorm(0, mean = 0, sd = sigma))
    }

    t_match <- regmatches(
        prior_str,
        regexec("^student_t\\(\\s*([^,]+)\\s*,\\s*([^,]+)\\s*,\\s*([^)]+)\\s*\\)$", prior_str)
    )[[1L]]
    if (length(t_match) == 4L) {
        df <- as.numeric(t_match[2L])
        mu <- as.numeric(t_match[3L])
        sigma <- as.numeric(t_match[4L])
        if (!isTRUE(mu == 0))
            cli::cli_abort(c(
                "Automatic {.arg parameter_prior} derivation requires a zero-centered slab prior for coefficient {.val {coef}}.",
                "i" = "Found {.val {prior_str}} (location = {mu}).",
                "i" = "Supply {.arg parameter_prior} explicitly."
            ))
        return(stats::dt(0, df = df) / sigma)
    }

    cli::cli_abort(c(
        "Cannot automatically derive {.arg parameter_prior} for coefficient {.val {coef}}.",
        "i" = "Found prior {.val {prior_str}}, which is not in the supported allowlist.",
        "i" = "Supported: {.val normal(0, s)} and {.val student_t(df, 0, s)}.",
        "i" = "Supply {.arg parameter_prior} explicitly."
    ))
}

# ──────────────────────────────────────────────────────────────────────────────
# Validation
# ──────────────────────────────────────────────────────────────────────────────

#' Validate sticky arguments for brm_pdmp
#'
#' @param sticky Logical.
#' @param can_stick Logical vector or NULL.
#' @param model_prior A bernoulli/betabernoulli object or NULL.
#' @param parameter_prior Numeric vector or NULL.
#' @param d Integer, total number of unconstrained parameters.
#' @param unc_names Character vector, unconstrained parameter names.
#' @param prior brms prior data frame.
#' @param subsampled Logical, whether subsampling is active.
#'
#' @return A list with validated `can_stick`, `model_prior`, and
#'   `parameter_prior` (all length `d`).
#'
#' @keywords internal
validate_brms_sticky <- function(sticky, can_stick, model_prior, parameter_prior,
                                  d, unc_names, prior, subsampled) {
    if (!isTRUE(sticky))
        return(list(sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL))

    if (subsampled)
        cli::cli_abort("Sticky variable selection with subsampled gradients is not yet supported. Use full-data gradients.")

    # model_prior is required
    if (is.null(model_prior) || !(is.bernoulli(model_prior) || is.betabernoulli(model_prior)))
        cli::cli_abort("{.arg model_prior} must be a {.cls bernoulli} or {.cls beta-bernoulli} object when {.arg sticky} is {.code TRUE}.")

    # Build can_stick
    can_stick_full <- map_can_stick(unc_names, user_can_stick = can_stick)
    n_stickable <- sum(can_stick_full)

    if (n_stickable == 0L)
        cli::cli_abort("No supported population-level coefficients found for variable selection. Check the model formula.")

    # Expand bernoulli prob
    if (is.bernoulli(model_prior)) {
        if (length(model_prior$prob) == 1L) {
            model_prior <- bernoulli(prob = rep(model_prior$prob, d))
        } else if (length(model_prior$prob) != d) {
            cli::cli_abort("For {.cls bernoulli} model_prior, {.arg prob} must be a scalar or length {d}.")
        }
    }

    # parameter_prior: derive automatically or validate user-supplied
    if (is.null(parameter_prior)) {
        parameter_prior <- derive_parameter_prior(prior, unc_names, can_stick_full)
    } else {
        if (length(parameter_prior) == n_stickable) {
            # Expand to full d-length vector
            pp_full <- rep(1.0, d)
            pp_full[which(can_stick_full)] <- parameter_prior
            parameter_prior <- pp_full
        } else if (length(parameter_prior) != d) {
            cli::cli_abort(c(
                "{.arg parameter_prior} must have length {n_stickable} (number of stickable coefficients) or {d} (total unconstrained parameters).",
                "x" = "Got length {length(parameter_prior)}."
            ))
        }
        validate_type(parameter_prior, type = "double", n = d, positive = TRUE)
    }

    list(
        sticky = TRUE,
        can_stick = can_stick_full,
        model_prior = model_prior,
        parameter_prior = parameter_prior
    )
}
