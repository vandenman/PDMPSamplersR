# Custom family definitions, link handling, codegen inspection, and Stan
# function generation for PDMP subsampling.

# -- Constants -----------------------------------------------------------------

.default_dpar_links <- list(
  sigma = "log",
  shape = "log"
)

.baked_in_links <- list(
  bernoulli   = c(logit = "logit"),
  poisson     = c(log = "log"),
  negbinomial = c(log = "log")
)

.family_configs <- list(
  bernoulli_logit = list(
    stan_dist     = "bernoulli_logit",
    dist_suffix   = "lpmf",
    response_type = "int",
    extra_dpars   = character(0)
  ),
  poisson_log = list(
    stan_dist     = "poisson_log",
    dist_suffix   = "lpmf",
    response_type = "int",
    extra_dpars   = character(0)
  ),
  gaussian = list(
    stan_dist     = "normal",
    dist_suffix   = "lpdf",
    response_type = "real",
    extra_dpars   = "sigma"
  ),
  negbinomial_log = list(
    stan_dist     = "neg_binomial_2_log",
    dist_suffix   = "lpmf",
    response_type = "int",
    extra_dpars   = "shape"
  )
)

# -- Family helpers ------------------------------------------------------------

subsampled_family_key <- function(fname, link) {
  baked <- .baked_in_links[[fname]]
  if (!is.null(baked) && !is.null(baked[[link]])) {
    paste0(fname, "_", baked[[link]])
  } else {
    fname
  }
}

subsampled_family_name <- function(fname, link) {
  baked <- .baked_in_links[[fname]]
  suffix <- if (!is.null(baked)) baked[[link]] else NULL
  if (!is.null(suffix)) {
    paste0("subsampled_", fname, "_", suffix)
  } else {
    paste0("subsampled_", fname)
  }
}

get_original_dpar_link <- function(family, dpar) {
  if (dpar == "mu") return(family$link)
  link <- .default_dpar_links[[dpar]]
  if (is.null(link))
    cli::cli_abort("Unknown default link for dpar {.val {dpar}} in family {.val {family$family}}.")
  link
}

get_subsampled_family_links <- function(family) {
  if (is.null(family$dpars))
    family <- brms::brmsfamily(family$family, link = family$link)
  dpars <- family$dpars
  links <- vapply(dpars, function(dp) get_original_dpar_link(family, dp), character(1))
  links[dpars == "mu"] <- "identity"
  unname(links)
}

.dpar_lower_bounds <- c(sigma = 0, shape = 0)

get_subsampled_family_bounds <- function(dpars) {
  vapply(dpars, function(dp) {
    b <- .dpar_lower_bounds[dp]
    if (is.na(b)) NA_real_ else unname(b)
  }, numeric(1), USE.NAMES = FALSE)
}

make_subsampled_family <- function(family) {
  if (is.character(family)) family <- do.call(family, list())
  if (is.null(family$dpars))
    family <- brms::brmsfamily(family$family, link = family$link)
  fname <- family$family
  link  <- family$link
  dpars <- family$dpars
  name  <- subsampled_family_name(fname, link)

  fam_key <- subsampled_family_key(fname, link)
  if (is.null(.family_configs[[fam_key]]))
    cli::cli_abort("Unsupported family for subsampling: {.val {fname}}(link = {.val {link}}).")

  type <- if (fname %in% c("bernoulli", "poisson", "negbinomial",
                            "binomial", "geometric")) "int" else "real"

  dpar_links <- get_subsampled_family_links(family)
  dpar_lb <- get_subsampled_family_bounds(dpars)

  brms::custom_family(
    name,
    dpars  = dpars,
    links  = dpar_links,
    lb     = dpar_lb,
    type   = type,
    vars   = "N",
    loop   = FALSE
  )
}

# -- Codegen inspection --------------------------------------------------------

inspect_subsampled_codegen <- function(formula, data, sub_family,
                                       prior = NULL, user_stanvars = NULL,
                                       sample_prior = "no", ...) {
  code <- brms::stancode(
    formula,
    data = data,
    family = sub_family,
    prior = prior,
    stanvars = user_stanvars,
    sample_prior = sample_prior,
    ...
  )

  model_block <- extract_named_block(code, "model")

  dpars <- sub_family$dpars
  non_mu_dpars <- setdiff(dpars, "mu")
  validate_subsampled_surface(model_block, non_mu_dpars)

  dpar_shapes <- list()
  for (dp in non_mu_dpars) {
    pattern <- paste0("vector\\[N\\]\\s+", dp, "\\s*=")
    if (grepl(pattern, model_block)) {
      dpar_shapes[[dp]] <- "vector"
    } else {
      dpar_shapes[[dp]] <- "scalar"
    }
  }

  list(dpar_shapes = dpar_shapes)
}

# -- Stan function generation --------------------------------------------------

make_subsampled_stan_functions <- function(family, codegen_info) {
  if (is.character(family)) family <- do.call(family, list())
  if (is.null(family$dpars))
    family <- brms::brmsfamily(family$family, link = family$link)
  fname <- family$family
  link  <- family$link

  fam_key <- subsampled_family_key(fname, link)
  config  <- .family_configs[[fam_key]]
  if (is.null(config))
    cli::cli_abort("Unsupported family for subsampling: {.val {fname}}(link = {.val {link}}).")

  func_name <- subsampled_family_name(fname, link)
  suffix    <- config$dist_suffix
  stan_dist <- config$stan_dist

  y_arg <- if (config$response_type == "int") "array[] int y" else "vector y"

  sig_parts  <- character(0)
  call_parts <- character(0)
  for (dp in config$extra_dpars) {
    shape <- codegen_info$dpar_shapes[[dp]] %||% "scalar"
    if (shape == "scalar") {
      sig_parts  <- c(sig_parts, paste0("real ", dp))
      call_parts <- c(call_parts, dp)
    } else {
      sig_parts  <- c(sig_parts, paste0("vector ", dp))
      call_parts <- c(call_parts, paste0(dp, "[i]"))
    }
  }

  all_sig  <- paste(c("vector mu", sig_parts, "int N_total"), collapse = ", ")
  all_call <- paste(c("mu[i]", call_parts), collapse = ", ")

  lines <- c(
    sprintf("  real %s_%s(%s, %s) {", func_name, suffix, y_arg, all_sig),
    "    int m = pdmp_get_subsample_size();",
    "    real ll = 0;",
    "    real scaling = 1.0 * N_total / m;",
    "    for (i in 1:m) {",
    "      int idx = pdmp_get_subsample_index(i);",
    sprintf("      ll += %s_%s(y[idx] | %s);", stan_dist, suffix, all_call),
    "    }",
    "    return ll * scaling;",
    "  }"
  )

  paste(lines, collapse = "\n")
}

# -- Main entry point: stanvars ------------------------------------------------

make_pdmp_stanvars <- function(formula, data, family, user_stanvars = NULL,
                                prior = NULL, sample_prior = "no", ...) {
  sub_family   <- make_subsampled_family(family)
  codegen_info <- inspect_subsampled_codegen(
    formula, data, sub_family,
    prior = prior, user_stanvars = user_stanvars,
    sample_prior = sample_prior, ...
  )

  sv_decls <- brms::stanvar(
    scode = paste(
      "  int pdmp_get_subsample_size();",
      "  int pdmp_get_subsample_index(int n);",
      "  matrix get_subsampled_Xc(matrix Xc_full);",
      "  array[] int get_subsampled_int_array(array[] int arr);",
      sep = "\n"
    ),
    block = "functions"
  )

  sv_fns <- brms::stanvar(
    scode = make_subsampled_stan_functions(family, codegen_info),
    block = "functions"
  )

  out <- sv_decls + sv_fns
  if (!is.null(user_stanvars)) out <- out + user_stanvars
  out
}

# -- Custom-family post-processing functions -----------------------------------

#' @export
log_lik_subsampled_bernoulli_logit <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  y  <- prep$data$Y[i]
  dbinom(y, size = 1, prob = plogis(mu), log = TRUE)
}

#' @export
posterior_predict_subsampled_bernoulli_logit <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  rbinom(length(mu), size = 1, prob = plogis(mu))
}

#' @export
posterior_epred_subsampled_bernoulli_logit <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  plogis(mu)
}

#' @export
log_lik_subsampled_poisson_log <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  y  <- prep$data$Y[i]
  dpois(y, lambda = exp(mu), log = TRUE)
}

#' @export
posterior_predict_subsampled_poisson_log <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  rpois(length(mu), lambda = exp(mu))
}

#' @export
posterior_epred_subsampled_poisson_log <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  exp(mu)
}

#' @export
log_lik_subsampled_gaussian <- function(i, prep) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y     <- prep$data$Y[i]
  dnorm(y, mean = mu, sd = sigma, log = TRUE)
}

#' @export
posterior_predict_subsampled_gaussian <- function(i, prep, ...) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  rnorm(length(mu), mean = mu, sd = sigma)
}

#' @export
posterior_epred_subsampled_gaussian <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  mu
}

#' @export
log_lik_subsampled_negbinomial_log <- function(i, prep) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  shape <- brms::get_dpar(prep, "shape", i = i)
  y     <- prep$data$Y[i]
  dnbinom(y, mu = exp(mu), size = shape, log = TRUE)
}

#' @export
posterior_predict_subsampled_negbinomial_log <- function(i, prep, ...) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  shape <- brms::get_dpar(prep, "shape", i = i)
  rnbinom(length(mu), mu = exp(mu), size = shape)
}

#' @export
posterior_epred_subsampled_negbinomial_log <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  exp(mu)
}
