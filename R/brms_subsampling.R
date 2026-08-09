
make_prior_standata <- function(sdata) {
  N_orig  <- sdata$N
  non_obs <- c("N", "K", "Kc", "prior_only", "means_X")
  prior   <- sdata
  prior$N <- 1L
  for (nm in setdiff(names(sdata), non_obs)) {
    val <- sdata[[nm]]
    if (is.matrix(val) && nrow(val) == N_orig) {
      prior[[nm]] <- val[1L, , drop = FALSE]
    } else if (!is.matrix(val) && is.atomic(val) && length(val) == N_orig) {
      prior[[nm]] <- array(val[1L])
    }
  }
  prior$prior_only <- 1L
  prior
}

make_opaque_deterministic_standata <- function(sdata) {
  prior <- sdata
  prior$prior_only <- 1L
  prior
}


subsampling_custom_code_reason <- function(stanvars) {
  if (is.null(stanvars)) return(NULL)
  model_stanvars <- Filter(
    function(x) identical(x$block %||% "", "model"), unclass(stanvars))
  text <- paste(vapply(model_stanvars, function(x) x$scode %||% "",
                       character(1)), collapse = "\n")
  if (grepl("\\bprior_only\\b", text, perl = TRUE)) {
    return("custom model code references prior_only and breaks deterministic closure")
  }
  assignment <- paste0(
    "\\b(?:mu[A-Za-z0-9_]*|Xc?(?:_[A-Za-z0-9_]+)?|offsets?)",
    "(?:\\s*\\[[^]]+\\])?\\s*(?:=(?!=)|\\+=|-=|\\*=|/=)"
  )
  if (grepl(assignment, text, perl = TRUE)) {
    return("custom model code modifies an accelerated predictor")
  }
  NULL
}

subsampling_eligibility <- function(formula, family, stanvars = NULL,
                                           sdata = NULL) {
  if (is.character(family)) family <- do.call(family, list())
  family_name <- family$family %||% ""
  link <- family$link %||% ""
  custom_reason <- subsampling_custom_code_reason(stanvars)
  if (!is.null(custom_reason)) {
    return(list(eligible = FALSE, reason = custom_reason))
  }
  supported_family <-
    (family_name %in% c("bernoulli", "binomial") && identical(link, "logit")) ||
    (identical(family_name, "poisson") && identical(link, "log")) ||
    (identical(family_name, "gaussian") && identical(link, "identity")) ||
    (family_name %in% c("categorical", "multinomial") && identical(link, "logit"))
  if (!supported_family) {
    return(list(eligible = FALSE,
                reason = "no certified subsampling likelihood provider exists for this family and link"))
  }

  terms <- tryCatch(
    brms::brmsterms(if (inherits(formula, c("brmsformula", "mvbrmsformula"))) formula else
      brms::bf(formula, family = family)),
    error = function(e) NULL
  )
  if (inherits(terms, "mvbrmsterms")) {
    response_families <- vapply(terms$terms,
      function(x) x$family$family %||% "", character(1))
    response_links <- vapply(terms$terms,
      function(x) x$family$link %||% "", character(1))
    affine <- !isTRUE(terms$rescor) &&
      all(response_families == "bernoulli") && all(response_links == "logit") &&
      all(vapply(terms$terms, function(response_terms) {
        mu <- response_terms$dpars$mu
        !is.null(mu$fe) && !isTRUE(mu$transform) &&
          length(response_terms$adforms %||% list()) == 0L &&
          length(setdiff(names(mu), c(
            "formula", "fe", "offset", "allvars", "family", "dpar",
            "resp", "respform", "adforms", "transform"
          ))) == 0L
      }, logical(1)))
    return(list(
      eligible = affine,
      reason = if (affine) NULL else
        "multivariate subsampling acceleration requires independent affine Bernoulli responses",
      has_intercept = TRUE,
      family = if (affine) "independent_bernoulli" else family_name,
      link = link,
      responses = terms$responses
    ))
  }
  if (is.null(terms) || isTRUE(terms$mv) || length(terms$nlpars) > 0L) {
    return(list(eligible = FALSE,
                reason = "the formula is not a univariate affine mean model"))
  }

  allowed_adforms <- switch(family_name,
    bernoulli = c("weights", "subset"),
    binomial = c("weights", "subset", "trials"),
    poisson = c("weights", "subset", "rate"),
    gaussian = c("weights", "subset", "se"),
    categorical = c("weights", "subset"),
    multinomial = c("weights", "subset", "trials"),
    character(0)
  )
  unsupported_adforms <- setdiff(names(terms$adforms %||% list()),
                                 allowed_adforms)
  if (length(unsupported_adforms) > 0L) {
    return(list(
      eligible = FALSE,
      reason = paste0(
        "no certified subsampling likelihood provider exists for response modifier: ",
        paste(unsupported_adforms, collapse = ", ")
      )
    ))
  }

  # Fail closed: the only likelihood geometry represented by the initial
  # envelope is the population-level model matrix plus an optional fixed
  # offset. Any other brms term component or addition-data modifier is outside
  # that affirmative description.
  allowed_top <- c(
    "formula", "family", "mv", "rescor", "mecor", "resp", "respform",
    "adforms", "dpars", "nlpars", "fdpars", "allvars"
  )
  extra_top <- setdiff(names(terms), allowed_top)
  if (length(extra_top) > 0L) {
    return(list(eligible = FALSE,
                reason = paste0("unsupported response metadata: ", paste(extra_top, collapse = ", "))))
  }

  # Categorical and multinomial category blocks are represented by X_mu*
  # matrices in standata rather than ordinary dpars metadata.
  if (family_name %in% c("categorical", "multinomial")) {
    if (length(terms$dpars) > 0L) {
      return(list(eligible = FALSE,
                  reason = "the categorical predictor geometry is not purely affine"))
    }
    return(list(eligible = TRUE, reason = NULL, has_intercept = TRUE,
                family = family_name, link = link))
  }

  required_dpars <- if (identical(family_name, "gaussian")) {
    intersect(c("mu", "sigma"), names(terms$dpars))
  } else "mu"
  if (!"mu" %in% required_dpars) {
    return(list(eligible = FALSE,
                reason = "the mean predictor is not represented by an affine block"))
  }
  for (dpar in required_dpars) {
    dpar_terms <- terms$dpars[[dpar]]
    allowed_dpar <- c(
      "formula", "fe", "offset", "allvars", "family", "dpar", "resp",
      "respform", "adforms", "transform"
    )
    extra_dpar <- setdiff(names(dpar_terms), allowed_dpar)
    transformed_ok <- identical(family_name, "gaussian") &&
      identical(dpar, "sigma") && identical(family$link_sigma %||% "log", "log")
    if (length(extra_dpar) > 0L || (isTRUE(dpar_terms$transform) && !transformed_ok) ||
        is.null(dpar_terms$fe)) {
      label <- if (length(extra_dpar)) paste(extra_dpar, collapse = ", ") else
        "transformed predictor"
      return(list(eligible = FALSE,
                  reason = paste0("unsupported predictor component: ", label)))
    }
  }
  mu_terms <- terms$dpars$mu
  allowed_mu <- c("formula", "fe", "offset", "allvars", "family", "dpar",
                  "resp", "respform", "adforms", "transform")
  extra_mu <- setdiff(names(mu_terms), allowed_mu)
  if (length(extra_mu) > 0L || isTRUE(mu_terms$transform) ||
      is.null(mu_terms$fe)) {
    label <- if (length(extra_mu)) paste(extra_mu, collapse = ", ") else "transformed predictor"
    return(list(eligible = FALSE,
                reason = paste0("unsupported predictor component: ", label)))
  }

  fe_terms <- stats::terms(mu_terms$fe)
  list(
    eligible = TRUE,
    reason = NULL,
    has_intercept = identical(attr(fe_terms, "intercept"), 1L),
    family = family_name,
    link = link
  )
}

subsampling_observation_multipliers <- function(sdata, family_name) {
  N <- as.integer(sdata$N)
  multiplier <- rep(1, N)
  if (!is.null(sdata$weights)) multiplier <- multiplier * as.numeric(sdata$weights)
  if (family_name %in% c("binomial", "multinomial")) {
    if (is.null(sdata$trials))
      cli::cli_abort("Binomial observation subsampling requires the brms trials vector.")
    multiplier <- multiplier * as.numeric(sdata$trials)
  }
  if (length(multiplier) != N || any(!is.finite(multiplier)) ||
      any(multiplier < 0)) {
    cli::cli_abort("Subsampling observation multipliers must be finite and nonnegative.")
  }
  multiplier
}

build_subsampling_compact_affine_design <- function(sdata, unc_names, has_intercept,
                                               dpar = "mu") {
  suffix <- if (identical(dpar, "mu")) "" else paste0("_", dpar)
  X_name <- paste0("X", suffix)
  X <- as.matrix(sdata[[X_name]])
  N <- as.integer(sdata$N)
  if (nrow(X) != N)
    cli::cli_abort("The brms population-level design matrix has the wrong number of rows.")

  columns <- list()
  indices <- integer()
  intercept_idx <- match(paste0("Intercept", suffix), unc_names)
  if (isTRUE(has_intercept)) {
    if (is.na(intercept_idx))
      cli::cli_abort("The brms intercept could not be mapped to a BridgeStan coordinate.")
    columns[[length(columns) + 1L]] <- rep(1, N)
    indices <- c(indices, intercept_idx)
    assign <- attr(X, "assign")
    if (ncol(X) < 1L || is.null(assign) || assign[[1L]] != 0L)
      cli::cli_abort("The brms model matrix does not contain the expected intercept column.")
    X_population <- X[, -1L, drop = FALSE]
  } else {
    if (!is.na(intercept_idx))
      cli::cli_abort("An unexpected BridgeStan intercept coordinate was found.")
    X_population <- X
  }
  if (ncol(X_population) > 0L) {
    X_affine <- if (isTRUE(has_intercept)) {
      means_name <- paste0("means_X", suffix)
      means <- if (!is.null(sdata[[means_name]]) && length(sdata[[means_name]]) == ncol(X_population)) {
        as.numeric(sdata[[means_name]])
      } else {
        colMeans(X_population)
      }
      sweep(X_population, 2L, means, "-")
    } else {
      X_population
    }
    b_idx <- match(paste0("b", suffix, ".", seq_len(ncol(X_affine))), unc_names)
    if (anyNA(b_idx)) {
      cli::cli_abort(c(
        "Could not map every population-level design column to a BridgeStan coordinate.",
        "i" = "Expected unconstrained names {.val {paste0('b.', seq_len(ncol(X_affine)))}}."
      ))
    }
    columns[[length(columns) + 1L]] <- X_affine
    indices <- c(indices, b_idx)
  }
  design <- if (length(columns)) do.call(cbind, columns) else matrix(0, N, 0L)
  order <- order(indices)
  list(
    design = design[, order, drop = FALSE],
    indices = as.integer(indices[order])
  )
}

build_subsampling_predictor_geometry <- function(sdata, unc_names, eligibility) {
  family <- eligibility$family
  N <- as.integer(sdata$N)
  d <- length(unc_names)
  dpars <- if (identical(family, "independent_bernoulli")) {
    eligibility$responses
  } else if (family %in% c("categorical", "multinomial")) {
    sub("^X_", "", grep("^X_mu", names(sdata), value = TRUE))
  } else if (identical(family, "gaussian") && !is.null(sdata$X_sigma)) {
    c("mu", "sigma")
  } else {
    "mu"
  }
  predictors <- lapply(dpars, function(dpar) {
    X_name <- if (identical(dpar, "mu")) "X" else paste0("X_", dpar)
    X <- sdata[[X_name]]
    has_intercept <- !is.null(attr(X, "assign")) &&
      length(attr(X, "assign")) > 0L && attr(X, "assign")[[1L]] == 0L
    build_subsampling_compact_affine_design(
      sdata, unc_names, has_intercept, dpar
    )
  })
  designs <- lapply(predictors, `[[`, "design")
  design_indices <- lapply(predictors, `[[`, "indices")

  offsets <- matrix(0, N, length(designs))
  for (k in seq_along(dpars)) {
    offset_name <- if (identical(dpars[[k]], "mu")) {
      "offsets"
    } else {
      paste0("offsets_", dpars[[k]])
    }
    value <- sdata[[offset_name]]
    if (!is.null(value)) offsets[, k] <- as.numeric(value)
  }
  if (identical(family, "poisson") && !is.null(sdata$denom)) {
    offsets[, 1L] <- offsets[, 1L] + log(as.numeric(sdata$denom))
  }

  if (identical(family, "gaussian") && length(designs) == 1L) {
    sigma_idx <- match("sigma", unc_names)
    if (!is.na(sigma_idx)) {
      designs[[2L]] <- matrix(1, N, 1L)
      design_indices[[2L]] <- as.integer(sigma_idx)
      offsets <- cbind(offsets, 0)
    }
  }

  response <- if (identical(family, "independent_bernoulli")) {
    do.call(cbind, lapply(dpars, function(resp) sdata[[paste0("Y_", resp)]]))
  } else as.matrix(sdata$Y)
  if (nrow(response) != N) response <- matrix(as.numeric(sdata$Y), N, 1L)
  se <- if (is.null(sdata$se)) numeric(0) else as.numeric(sdata$se)
  list(
    designs = designs,
    design_indices = design_indices,
    dimension = d,
    offsets = offsets,
    response = response,
    se = se
  )
}

#' Show the Stan code used by brm_pdmp
#'
#' This is a thin wrapper around [brms::stancode()] using the same model
#' arguments accepted by [brm_pdmp()]. Observation subsampling uses the ordinary
#' Stan model plus analytic Julia-side family providers, so no rewritten Stan
#' variant is generated.
#'
#' @inheritParams brm_pdmp
#' @return A character string containing the generated Stan program.
#' @export
brm_stancode <- function(
    formula, data, family = gaussian(), prior = NULL,
    stanvars = NULL, sample_prior = "no", ...
) {
  if (!requireNamespace("brms", quietly = TRUE))
    cli::cli_abort("Package {.pkg brms} is required for {.fn brm_stancode}.")
  brms::stancode(formula, data = data, family = family, prior = prior,
                  stanvars = stanvars, sample_prior = sample_prior, ...)
}

#' Show the Stan data used by brm_pdmp
#'
#' This is a thin wrapper around [brms::standata()]. Observation subsampling derives
#' its deterministic and family geometry from this ordinary standata object.
#'
#' @inheritParams brm_pdmp
#' @return A named list containing the generated Stan data.
#' @export
brm_standata <- function(
    formula, data, family = gaussian(), prior = NULL,
    stanvars = NULL, sample_prior = "no", ...
) {
  if (!requireNamespace("brms", quietly = TRUE))
    cli::cli_abort("Package {.pkg brms} is required for {.fn brm_standata}.")
  brms::standata(formula, data = data, family = family, prior = prior,
                 stanvars = stanvars, sample_prior = sample_prior, ...)
}
