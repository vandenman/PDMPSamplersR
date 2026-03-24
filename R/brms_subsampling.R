subset_standata <- function(sdata, indices) {
  N_orig  <- sdata$N
  non_obs <- c("N", "K", "Kc", "prior_only", "means_X")
  sub     <- sdata
  sub$N   <- length(indices)
  for (nm in setdiff(names(sdata), non_obs)) {
    val <- sdata[[nm]]
    if (is.matrix(val) && nrow(val) == N_orig) {
      sub[[nm]] <- val[indices, , drop = FALSE]
    } else if (!is.matrix(val) && is.atomic(val) && length(val) == N_orig) {
      sub_val <- val[indices]
      if (length(sub_val) == 1L) sub_val <- array(sub_val)
      sub[[nm]] <- sub_val
    }
  }
  sub
}

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

inject_ext_cpp_stancode <- function(stancode) {
  y_is_int <- grepl("array\\[N\\] int Y;", stancode)

  if (y_is_int) {
    func_decls <- paste(
      "  array[] int get_subsampled_Y_int(array[] int Y_full);",
      "  matrix get_subsampled_Xc(matrix Xc_full);",
      sep = "\n"
    )
    y_getter <- "get_subsampled_Y_int"
  } else {
    func_decls <- paste(
      "  vector get_subsampled_Y_real(vector Y_full);",
      "  matrix get_subsampled_Xc(matrix Xc_full);",
      sep = "\n"
    )
    y_getter <- "get_subsampled_Y_real"
  }

  stancode <- sub(
    "(functions\\s*\\{)\n",
    paste0("\\1\n", func_decls, "\n"),
    stancode
  )

  gsub(
    "(_glm_lp[dm]f\\()([^|]+)(\\|\\s*)Xc,",
    paste0("\\1", y_getter, "(\\2) \\3get_subsampled_Xc(Xc),"),
    stancode
  )
}

hpp_path <- function() {
  system.file("stan", "pdmp_subsample.hpp", package = "PDMPSamplersR")
}

brms_has_subsample <- function() {
  "subsample" %in% names(formals(brms:::stancode.default))
}

pdmp_subsampling <- function(stancode) {
  y_is_int <- grepl("array\\[N\\] int Y;", stancode)
  y_getter <- if (y_is_int) "get_subsampled_Y_int" else "get_subsampled_Y_real"
  brms::subsampling(
    size_fn = "pdmp_get_subsample_size",
    index_fn = "pdmp_get_subsample_index",
    wrap = list(Y = y_getter, Xc = "get_subsampled_Xc")
  )
}

pdmp_function_decls <- function(stancode) {
  y_is_int <- grepl("array\\[N\\] int Y;", stancode)
  y_decl <- if (y_is_int) {
    "  array[] int get_subsampled_Y_int(array[] int Y_full);"
  } else {
    "  vector get_subsampled_Y_real(vector Y_full);"
  }
  paste(
    "  int pdmp_get_subsample_size();",
    "  int pdmp_get_subsample_index(int n);",
    y_decl,
    "  matrix get_subsampled_Xc(matrix Xc_full);",
    sep = "\n"
  )
}

make_ext_cpp_stancode <- function(standard_code, formula, data, family,
                                  prior, stanvars, sample_prior, ...) {
  if (brms_has_subsample()) {
    sub <- pdmp_subsampling(standard_code)
    decls <- pdmp_function_decls(standard_code)
    sv <- brms::stanvar(scode = decls, block = "functions")
    stanvars <- if (is.null(stanvars)) sv else stanvars + sv
    brms::stancode(formula, data = data, family = family,
                   prior = prior, stanvars = stanvars,
                   sample_prior = sample_prior,
                   subsample = sub, ...)
  } else {
    inject_ext_cpp_stancode(standard_code)
  }
}

#' Show the Stan code used by brm_pdmp
#'
#' Returns the Stan code that [brm_pdmp()] would compile, after any
#' subsampling modifications. Without `subsample_size`, this is
#' identical to `brms::stancode()`. With `subsample_size`, two
#' variants are returned:
#' \describe{
#'   \item{standard}{The standard model (used for the full-data and
#'     prior-only models).}
#'   \item{ext_cpp}{The model with likelihood rewritten to use
#'     external C++ subsetting functions (used for the
#'     subsampled-gradient model).}
#' }
#'
#' @inheritParams brm_pdmp
#'
#' @return A character string when `subsample_size` is NULL, or a
#'   named list with elements `standard` and `ext_cpp`.
#' @export
brm_stancode <- function(
    formula, data, family = gaussian(),
    prior = NULL,
    subsample_size = NULL,
    stanvars = NULL, sample_prior = "no",
    ...
) {
  if (!requireNamespace("brms", quietly = TRUE))
    cli::cli_abort("Package {.pkg brms} is required for {.fn brm_stancode}.")

  scode <- brms::stancode(formula, data = data, family = family,
                          prior = prior, stanvars = stanvars,
                          sample_prior = sample_prior, ...)
  if (is.null(subsample_size)) return(scode)

  list(
    standard = scode,
    ext_cpp  = make_ext_cpp_stancode(scode, formula, data, family,
                                     prior, stanvars, sample_prior, ...)
  )
}

#' Show the Stan data used by brm_pdmp
#'
#' Returns the Stan data that [brm_pdmp()] would pass to BridgeStan,
#' after any subsampling modifications. Without `subsample_size`, this
#' is identical to `brms::standata()`. With `subsample_size`, returns
#' a named list with three elements:
#' \describe{
#'   \item{full}{Full dataset.}
#'   \item{prior}{Prior-only dummy data (`N=1`, `prior_only=1`).}
#'   \item{subsample}{An initial subsample of `subsample_size` rows.}
#' }
#'
#' @inheritParams brm_pdmp
#' @param indices Optional integer vector of observation indices for
#'   the initial subsample. If NULL (default), a random sample of
#'   size `subsample_size` is drawn.
#'
#' @return A list (Stan data) when `subsample_size` is NULL, or a
#'   named list with elements `full`, `prior`, and `subsample`.
#' @export
brm_standata <- function(
    formula, data, family = gaussian(),
    prior = NULL,
    subsample_size = NULL,
    indices = NULL,
    stanvars = NULL, sample_prior = "no",
    ...
) {
  if (!requireNamespace("brms", quietly = TRUE))
    cli::cli_abort("Package {.pkg brms} is required for {.fn brm_standata}.")

  sdata <- brms::standata(formula, data = data, family = family,
                          prior = prior, stanvars = stanvars,
                          sample_prior = sample_prior, ...)

  if (is.null(subsample_size)) return(sdata)

  N <- nrow(data)
  subsample_size <- as.integer(subsample_size)
  if (subsample_size >= N)
    cli::cli_abort("{.arg subsample_size} ({subsample_size}) must be less than {.code nrow(data)} ({N}).")

  sdata_prior <- make_prior_standata(sdata)

  if (is.null(indices)) {
    indices <- sample.int(N, subsample_size)
  } else {
    if (length(indices) != subsample_size)
      cli::cli_abort("{.arg indices} length ({length(indices)}) must equal {.arg subsample_size} ({subsample_size}).")
  }
  sdata_sub <- subset_standata(sdata, indices)

  list(full = sdata, prior = sdata_prior, subsample = sdata_sub)
}
