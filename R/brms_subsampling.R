fix_brms_stancode <- function(stancode) {
  lines <- strsplit(stancode, "\n")[[1]]
  n <- length(lines)

  td_start <- grep("^transformed data\\s*\\{", lines)
  if (length(td_start) == 0) return(stancode)

  means_decl <- grep("^\\s*vector\\[Kc\\] means_X;", lines)
  means_assign <- grep("means_X\\[i - 1\\] = mean\\(X\\[, i\\]\\);", lines)
  if (length(means_decl) == 0 || length(means_assign) == 0) return(stancode)

  data_close <- grep("^\\s*int prior_only;", lines)
  if (length(data_close) == 0) {
    data_close <- grep("^\\}", lines)
    data_close <- data_close[data_close < td_start[1]]
    if (length(data_close) == 0) return(stancode)
    data_close <- max(data_close)
  } else {
    data_close <- data_close[1]
  }

  new_lines <- character(0)
  skip <- integer(0)
  skip <- c(means_decl[1], means_assign[1])
  inserted <- FALSE

  for (i in seq_len(n)) {
    if (i %in% skip) next
    new_lines <- c(new_lines, lines[i])
    if (i == data_close && !inserted) {
      new_lines <- c(new_lines, "  vector[Kc] means_X;  // column means of X (fixed for subsampling)")
      inserted <- TRUE
    }
  }

  paste(new_lines, collapse = "\n")
}

subset_standata <- function(sdata, indices, means_X) {
  sub <- sdata
  sub$N <- length(indices)
  Y_sub <- sdata$Y[indices]
  if (length(Y_sub) == 1L) Y_sub <- array(Y_sub)
  sub$Y <- Y_sub
  sub$X <- sdata$X[indices, , drop = FALSE]
  sub$means_X <- means_X
  sub
}

make_prior_standata <- function(sdata, means_X) {
  K <- sdata$K
  prior <- list(
    N = 1L,
    Y = array(if (is.integer(sdata$Y)) 0L else 0.0),
    K = K,
    Kc = sdata$Kc,
    X = matrix(0, nrow = 1, ncol = K),
    means_X = means_X,
    prior_only = 1L
  )
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

  stancode <- gsub(
    "(_glm_lp[dm]f\\()([^|]+)(\\|\\s*)Xc,",
    paste0("\\1", y_getter, "(\\2) \\3get_subsampled_Xc(Xc),"),
    stancode
  )
  stancode
}

hpp_path <- function() {
  system.file("stan", "pdmp_subsample.hpp", package = "PDMPSamplersR")
}

#' Show the Stan code used by brm_pdmp
#'
#' Returns the Stan code that [brm_pdmp()] would compile, after any
#' subsampling modifications. Without `subsample_size`, this is
#' identical to `brms::stancode()`. With `subsample_size`, two
#' variants are returned:
#' \describe{
#'   \item{standard}{The standard model with `means_X` moved to data
#'     (used for the full-data and prior-only models).}
#'   \item{ext_cpp}{The model with external C++ subsetting functions
#'     injected (used for the subsampled-gradient model).}
#' }
#'
#' @inheritParams brm_pdmp
#'
#' @return A character string when `subsample_size` is NULL, or a
#'   named list with elements `standard` and `ext_cpp`.
#' @export
pdmp_stancode <- function(
    formula, data, family = gaussian(),
    prior = NULL,
    subsample_size = NULL,
    stanvars = NULL, sample_prior = "no",
    ...
) {
  if (!requireNamespace("brms", quietly = TRUE))
    cli::cli_abort("Package {.pkg brms} is required for {.fn pdmp_stancode}.")

  scode <- brms::stancode(formula, data = data, family = family,
                          prior = prior, stanvars = stanvars,
                          sample_prior = sample_prior, ...)
  if (is.null(subsample_size)) return(scode)

  list(
    standard = fix_brms_stancode(scode),
    ext_cpp  = inject_ext_cpp_stancode(fix_brms_stancode(scode))
  )
}

#' Show the Stan data used by brm_pdmp
#'
#' Returns the Stan data that [brm_pdmp()] would pass to BridgeStan,
#' after any subsampling modifications. Without `subsample_size`, this
#' is identical to `brms::standata()`. With `subsample_size`, returns
#' a named list with three elements:
#' \describe{
#'   \item{full}{Full dataset with `means_X` added.}
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
pdmp_standata <- function(
    formula, data, family = gaussian(),
    prior = NULL,
    subsample_size = NULL,
    indices = NULL,
    stanvars = NULL, sample_prior = "no",
    ...
) {
  if (!requireNamespace("brms", quietly = TRUE))
    cli::cli_abort("Package {.pkg brms} is required for {.fn pdmp_standata}.")

  sdata <- brms::standata(formula, data = data, family = family,
                          prior = prior, stanvars = stanvars,
                          sample_prior = sample_prior, ...)

  if (is.null(subsample_size)) return(sdata)

  N <- nrow(data)
  subsample_size <- as.integer(subsample_size)
  if (subsample_size >= N)
    cli::cli_abort("{.arg subsample_size} ({subsample_size}) must be less than {.code nrow(data)} ({N}).")

  if (any(grepl("^(J|Z)_", names(sdata))))
    cli::cli_abort(c(
      "Subsampled gradients currently support fixed-effects models only.",
      "i" = "Remove random effects from the formula, or omit {.arg subsample_size} to use full-data gradients."
    ))

  means_X <- array(colMeans(sdata$X[, -1, drop = FALSE]))
  sdata_full <- sdata
  sdata_full$means_X <- means_X
  sdata_prior <- make_prior_standata(sdata, means_X)

  if (is.null(indices)) {
    indices <- sample.int(N, subsample_size)
  } else {
    if (length(indices) != subsample_size)
      cli::cli_abort("{.arg indices} length ({length(indices)}) must equal {.arg subsample_size} ({subsample_size}).")
  }
  sdata_sub <- subset_standata(sdata, indices, means_X)

  list(full = sdata_full, prior = sdata_prior, subsample = sdata_sub)
}
