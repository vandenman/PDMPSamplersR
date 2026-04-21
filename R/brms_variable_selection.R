# Sticky variable selection helpers for brm_pdmp()
#
# Maps brms formula terms to unconstrained BridgeStan coordinates and
# derives per-coefficient kappa values from brms prior specifications.

# ──────────────────────────────────────────────────────────────────────────────
# can_stick mapping
# ──────────────────────────────────────────────────────────────────────────────

supported_b_coef_names <- function(fe_names, formula = NULL, data = NULL) {
    if (!is.character(fe_names))
        cli::cli_abort("{.arg fe_names} must be a character vector.")
    if (length(fe_names) == 0L)
        return(character(0))
    if (is.null(formula))
        return(paste0("b.", fe_names))
    if (is.null(data))
        cli::cli_abort("{.arg data} must be provided when {.arg formula} is supplied.")

    main_formula <- .extract_main_formula(formula)
    formula_text <- paste(deparse(formula), collapse = " ")
    tilde_pos <- gregexpr("~", formula_text, fixed = TRUE)[[1L]]
    n_tilde <- if (length(tilde_pos) == 1L && tilde_pos[1L] == -1L) 0L else length(tilde_pos)
    if (n_tilde > 1L)
        cli::cli_abort(c(
            "Formula-based sticky auto-mapping supports only a single-response formula.",
            "i" = "Distributional, multi-response, and non-linear formula extensions are not supported."
        ))

    term_labels <- attr(stats::terms(main_formula), "term.labels")
    if (length(term_labels) == 0L)
        return(character(0))

    unsupported_terms <- term_labels[grepl("[:|]", term_labels) | grepl("\\(|\\)|\\[|\\]", term_labels)]
    if (length(unsupported_terms) > 0L)
        cli::cli_abort(c(
            "Formula-based sticky auto-mapping supports only simple numeric main effects.",
            "x" = "Unsupported terms: {.val {unsupported_terms}}.",
            "i" = "Interactions, splines, random effects, and transformed terms are not supported."
        ))

    unknown_terms <- setdiff(term_labels, names(data))
    if (length(unknown_terms) > 0L)
        cli::cli_abort(c(
            "Formula-based sticky auto-mapping requires direct data-column predictors.",
            "x" = "Terms not found as plain data columns: {.val {unknown_terms}}.",
            "i" = "Use simple numeric predictors like {.code y ~ x1 + x2}."
        ))

    non_numeric_terms <- term_labels[!vapply(term_labels, function(lbl) is.numeric(data[[lbl]]), logical(1L))]
    if (length(non_numeric_terms) > 0L)
        cli::cli_abort(c(
            "Formula-based sticky auto-mapping supports numeric predictors only.",
            "x" = "Non-numeric terms: {.val {non_numeric_terms}}.",
            "i" = "Factors and grouped structures are not supported."
        ))

    if (!setequal(fe_names, term_labels))
        cli::cli_abort(c(
            "Formula-based sticky auto-mapping requires one fixed-effect coefficient per predictor term.",
            "x" = "Formula terms: {.val {term_labels}}.",
            "x" = "Fixed-effect coefficients: {.val {fe_names}}.",
            "i" = "This mismatch indicates an unsupported model structure for sticky auto-mapping."
        ))

    paste0("b.", fe_names)
}

.extract_main_formula <- function(formula) {
    if (inherits(formula, "formula"))
        return(formula)
    if (is.list(formula) && !is.null(formula$formula) && inherits(formula$formula, "formula"))
        return(formula$formula)

    cli::cli_abort(c(
        "Unsupported {.arg formula} type for formula-based sticky auto-mapping.",
        "i" = "Use a simple single-response formula like {.code y ~ x1 + x2}."
    ))
}

#' Build logical can_stick vector from brms metadata
#'
#' Determines which unconstrained BridgeStan coordinates correspond to
#' supported population-level coefficients (excluding intercept, scale,
#' and shape parameters).
#'
#' @param unc_names Character vector from `BridgeStan::param_unc_names()`.
#' @param supported_coef_names Character vector of supported coefficient names
#'   on the unconstrained scale (e.g. `b.x1`) derived from brms metadata.
#' @param user_can_stick Optional logical vector (same length as supported
#'   non-intercept `b` coefficients) to override the default all-TRUE.
#'
#' @return A logical vector of length `length(unc_names)`.
#'
#' @keywords internal
map_can_stick <- function(unc_names, supported_coef_names = NULL, user_can_stick = NULL) {
    stickable_names <- stickable_coef_names(unc_names, supported_coef_names = supported_coef_names)
    stickable_indices <- match(stickable_names, unc_names)
    out <- rep(FALSE, length(unc_names))
    if (length(stickable_indices) == 0L)
        return(out)

    if (!is.null(user_can_stick)) {
        if (!is.logical(user_can_stick))
            cli::cli_abort("{.arg can_stick} must be a logical vector.")

        if (!is.null(names(user_can_stick))) {
            if (any(is.na(names(user_can_stick)) | !nzchar(names(user_can_stick))))
                cli::cli_abort("Named {.arg can_stick} vectors cannot contain empty or NA names.")

            # Named vector: match by name, accepting names with or without b. prefix
            stickable_short <- sub("^b\\.", "", stickable_names)
            norm_user_names <- sub("^b\\.", "", names(user_can_stick))
            if (anyDuplicated(norm_user_names))
                cli::cli_abort("Named {.arg can_stick} cannot contain duplicate coefficient names.")

            unmatched <- setdiff(norm_user_names, stickable_short)
            if (length(unmatched) > 0L)
                cli::cli_abort(c(
                    "{.arg can_stick} contains names not found among stickable coefficients.",
                    "x" = "Unknown names: {.val {unmatched}}.",
                    "i" = "Stickable coefficients: {.val {stickable_short}}."
                ))

            match_idx <- match(norm_user_names, stickable_short)
            out[stickable_indices[match_idx]] <- user_can_stick
        } else {
            if (length(user_can_stick) != length(stickable_indices))
                cli::cli_abort(c(
                    "{.arg can_stick} must have length {length(stickable_indices)} (the number of supported non-intercept population-level coefficients).",
                    "i" = "Stickable coefficients: {.val {stickable_names}}.",
                    "x" = "Got length {length(user_can_stick)}."
                ))
            out[stickable_indices] <- user_can_stick
        }
    } else {
        out[stickable_indices] <- TRUE
    }

    out
}

#' Return the names of stickable coefficients
#'
#' @param unc_names Character vector from `BridgeStan::param_unc_names()`.
#' @param supported_coef_names Optional character vector of supported
#'   coefficient names from brms metadata. If omitted, uses legacy raw
#'   name matching over `unc_names`.
#' @return Character vector of stickable unconstrained parameter names.
#' @keywords internal
stickable_coef_names <- function(unc_names, supported_coef_names = NULL) {
    if (is.null(supported_coef_names)) {
        b_indices <- grep("^b\\.", unc_names)
        intercept_indices <- grep("^b\\..*[Ii]ntercept", unc_names)
        stickable_indices <- setdiff(b_indices, intercept_indices)
        return(unc_names[stickable_indices])
    }

    if (!is.character(supported_coef_names))
        cli::cli_abort("{.arg supported_coef_names} must be a character vector.")
    if (anyDuplicated(supported_coef_names))
        cli::cli_abort("{.arg supported_coef_names} cannot contain duplicates.")
    if (length(supported_coef_names) == 0L)
        return(character(0))

    idx <- match(supported_coef_names, unc_names)
    missing <- supported_coef_names[is.na(idx)]
    if (length(missing) > 0L)
        cli::cli_abort(c(
            "Could not align supported brms coefficients to unconstrained parameter names.",
            "x" = "Missing from {.arg unc_names}: {.val {missing}}.",
            "i" = "This likely indicates drift between brms naming and BridgeStan parameter naming."
        ))

    unc_names[idx]
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
#' @param supported_coef_names Character vector of supported coefficient names
#'   derived from brms metadata.
#' @param prior brms prior data frame.
#' @param subsampled Logical, whether subsampling is active.
#'
#' @return A list with validated `can_stick`, `model_prior`, and
#'   `parameter_prior` (all length `d`).
#'
#' @keywords internal
validate_brms_sticky <- function(sticky, can_stick, model_prior, parameter_prior,
                                  d, unc_names, supported_coef_names = NULL, prior, subsampled) {
    if (!isTRUE(sticky))
        return(list(sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL))

    if (isTRUE(subsampled))
        cli::cli_warn(c(
            "Sticky dynamics with subsampled gradients are experimental and do not yet provide a runtime advantage.",
            "i" = "Benchmarking shows the subsampled path is ~12-16x slower than full-data sticky sampling with default settings.",
            "i" = "Further investigation is planned at the Julia package level. Use with caution."
        ))

    # model_prior is required
    if (is.null(model_prior) || !(is.bernoulli(model_prior) || is.betabernoulli(model_prior)))
        cli::cli_abort("{.arg model_prior} must be a {.cls bernoulli} or {.cls beta-bernoulli} object when {.arg sticky} is {.code TRUE}.")

    supported_unc_names <- stickable_coef_names(unc_names, supported_coef_names = supported_coef_names)
    if (length(supported_unc_names) == 0L)
        cli::cli_abort("No supported population-level coefficients found for variable selection. Check the model formula.")

    # Build can_stick
    can_stick_full <- map_can_stick(
        unc_names,
        supported_coef_names = supported_coef_names,
        user_can_stick = can_stick
    )
    n_stickable <- sum(can_stick_full)

    if (n_stickable == 0L)
        cli::cli_abort(c(
            "No coefficients are currently selected for sticky variable selection.",
            "i" = "Supported coefficients: {.val {supported_unc_names}}.",
            "i" = "Check {.arg can_stick}."
        ))

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
        supported_coef_names = supported_unc_names,
        model_prior = model_prior,
        parameter_prior = parameter_prior
    )
}
