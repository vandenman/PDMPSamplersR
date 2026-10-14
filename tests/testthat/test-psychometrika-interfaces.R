test_that("public subsampling API has no legacy terminology", {
  legacy_term <- paste0("mark", "ed")
  exports <- getNamespaceExports("PDMPSamplersR")
  expect_false(any(grepl(legacy_term, exports, ignore.case = TRUE)))
  exported_functions <- Filter(function(name) {
    is.function(getExportedValue("PDMPSamplersR", name))
  }, exports)
  public_arguments <- unlist(lapply(exported_functions, function(name) {
    names(formals(getExportedValue("PDMPSamplersR", name)))
  }), use.names = FALSE)
  expect_false(any(grepl(legacy_term, public_arguments, ignore.case = TRUE)))

  spec <- stan_subsampling(
    3L, 1L, list(prior_only = 1L), stan_residual_envelope(rep(1, 3L))
  )
  expect_identical(class(spec), "stan_subsampling")
  expect_false(any(grepl(legacy_term, c(class(spec), names(spec)),
                         ignore.case = TRUE)))
})

test_that("subsampling Stan specification validates certified envelopes", {
  envelope <- stan_residual_envelope(matrix(c(1, 2, 3), nrow = 1))
  subsampling <- stan_subsampling(
    n_observations = 3,
    subsample_size = 1,
    prior_standata = list(prior_only = 1L),
    residual_envelope = envelope
  )
  expect_s3_class(subsampling, "stan_subsampling")
  expect_equal(subsampling$subsample_size, 1L)
  full_batch <- stan_subsampling(3, 3, list(), envelope)
  expect_equal(full_batch$subsample_size, 3L)
  expect_error(
    stan_subsampling(3, 1, list(), list(weights = matrix(1, 1, 3))),
    "certified"
  )
})

test_that("unsupported custom-Stan controls fail before model construction", {
  subsampling <- stan_subsampling(
    4L, 2L, testthat::test_path("..", "stan", "subsampling_gaussian_prior.json"),
    stan_residual_envelope(rep(1, 4L))
  )
  base <- list(
    path_to_stanmodel = testthat::test_path(
      "..", "stan", "subsampling_gaussian.stan"),
    standata = testthat::test_path(
      "..", "stan", "subsampling_gaussian_full.json"),
    subsampling = subsampling, T = 0.01, show_progress = FALSE
  )
  expect_error(do.call(pdmp_sample_from_stanmodel,
                       c(base, list(curvature_backend = "exact"))),
               "no deterministic exact HVP")
  expect_error(do.call(pdmp_sample_from_stanmodel, c(base, list(
    support_boundary = support_boundary_control(mode = "line_search")
  ))), "requires.*support_boundary.*error")
})

test_that("OMRF provider specifies the whole-person direct HCV envelope", {
  X <- rbind(c(0, 1, 2), c(2, 0, 1), c(1, 2, 0))
  envelope <- omrf_residual_envelope(
    X, seen = c(3, 3, 3),
    thresholds = "thresholds_0",
    interactions = "interactions_0"
  )
  expect_s3_class(envelope, "omrf_residual_envelope")
  expect_equal(ncol(envelope$weights), nrow(X))
  expect_equal(dim(envelope$node_weights), c(ncol(X), nrow(X)))
  expect_equal(unname(envelope$edge_order),
               rbind(c(1L, 2L), c(1L, 3L), c(2L, 3L)))
  expect_error(
    omrf_residual_envelope(X, c(2, 3, 3), "a", "b"),
    "encoded"
  )
  expect_error(
    omrf_residual_envelope(X, c(3, 3, 3), "a", "b", geometry = "pattern_local"),
    "unused argument"
  )
  expect_false(omrf_residual_envelope(
    X, c(3, 3, 3), "a", "b",
    exact_selected_refinement = FALSE)$exact_selected_refinement)
  expect_error(
    omrf_residual_envelope(X, c(3, 3, 3), "a", "b",
                           signed_only_person_roof = TRUE),
    "person_signed_roof"
  )
  expect_error(
    omrf_residual_envelope(X, c(3, 3, 3), "a", "b",
                           shared_current_categorical = TRUE),
    "joint groups"
  )
})

test_that("OMRF analytic prior routing prefers an explicit independent prior", {
  X <- matrix(0L, nrow = 2L, ncol = 2L)
  spec <- omrf_residual_envelope(
    X, seen = c(2L, 2L),
    thresholds = "thresholds_0", interactions = "interactions_0")
  standata <- list(
    prior_threshold_alpha = 1,
    prior_threshold_beta = 1,
    prior_cauchy_scale = 1,
    # These fields are also present in the shared OMRF simulation data, but
    # belong to the dependent node-scale Stan model.
    slab_base_scale = 1,
    prior_global_scale_sd = 1,
    prior_node_scale_sd = 1
  )
  attach_prior <- get(".attach_omrf_analytic_prior", asNamespace("PDMPSamplersR"))
  resolved <- attach_prior(spec, standata)
  expect_identical(resolved$analytic_prior$type, "cauchy")
  expect_null(resolved$analytic_prior$global_scale)
  expect_null(resolved$analytic_prior$node_scale)
})

test_that("parameter blocks resolve to ordered unconstrained coordinates", {
  names <- c("thresholds_0.1", "thresholds_0.2",
             "interactions_0.1", "log_tau")
  expect_equal(
    PDMPSamplersR:::.resolve_unconstrained_spec(
      c("thresholds_0", "log_tau"), names),
    c(1L, 2L, 4L)
  )
  expect_error(
    PDMPSamplersR:::.resolve_unconstrained_spec(
      c("thresholds_0", "thresholds_0.1"), names),
    "duplicates"
  )
})

test_that("custom-Stan anchor-bank controls validate before backend setup", {
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", list(), warmup_adaptation_interval = -0.1),
    "finite nonnegative"
  )
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", list(), warmup_adaptation_interval = Inf),
    "finite nonnegative"
  )
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", list(), subsampling_anchor_updates = 1L,
      t_warmup = 1),
    "require.*subsampling"
  )
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", list(), subsampling_use_anchor_bank = TRUE),
    "requires positive"
  )
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", list(), subsampling_anchor_bank_capacity = 0L),
    "positive integer"
  )
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", list(),
      subsampling_main_anchor_refresh_distance = 0),
    "positive or Inf"
  )
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", list(),
      subsampling_main_anchor_refresh_distance = 0.5),
    "requires.*anchor_bank"
  )
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", list(),
      subsampling_require_no_main_full_data = NA),
    "must be TRUE or FALSE"
  )
})

test_that("log-linear Gaussian scale slab validates compact design", {
  slab <- loglinear_gaussian_scale_slab(
    log_base_sd = log(0.5),
    logscale = c("log_tau", "log_lambda"),
    logscale_design = matrix(c(1, 0.5, 1, 0.5), nrow = 2),
    coef = "interactions_0"
  )
  expect_s3_class(slab, "dependent_slab_prior")
  expect_equal(slab$type, "loglinear_gaussian_scale")
  expect_error(
    loglinear_gaussian_scale_slab(0, 1L, matrix(1, 2, 2)),
    "Columns"
  )
})

test_that("log-linear Gaussian scale slab canonicalizes supported Matrix inputs", {
  skip_if_not_installed("Matrix")
  csc <- Matrix::sparseMatrix(
    i = c(1L, 2L, 2L), j = c(1L, 1L, 2L), x = c(1, 0.5, 1),
    dims = c(2L, 2L)
  )
  fixtures <- list(
    csc = csc,
    triplet = Matrix::sparseMatrix(
      i = c(1L, 2L, 2L), j = c(1L, 1L, 2L), x = c(1, 0.5, 1),
      dims = c(2L, 2L), repr = "T"),
    row_sparse = Matrix::sparseMatrix(
      i = c(1L, 2L, 2L), j = c(1L, 1L, 2L), x = c(1, 0.5, 1),
      dims = c(2L, 2L), repr = "R"),
    symmetric = Matrix::forceSymmetric(Matrix::sparseMatrix(
      i = c(1L, 1L, 2L), j = c(1L, 2L, 2L), x = c(1, 0.5, 1),
      dims = c(2L, 2L)))
  )
  for (design in fixtures) {
    slab <- loglinear_gaussian_scale_slab(
      0, logscale = c("global", "node"),
      logscale_design = design, coef = "interactions_0"
    )
    transport <- slab$logscale_design
    expect_identical(transport$storage, "sparse_rows")
    expect_identical(transport$dims, c(2L, 2L))
    reconstructed <- matrix(0, 2L, 2L)
    reconstructed[cbind(transport$i, transport$j)] <- transport$x
    expect_equal(reconstructed, as.matrix(design))
  }

  dense <- loglinear_gaussian_scale_slab(
    0, c("global", "node"), matrix(c(1, 0.5, 0, 1), 2L))
  expect_identical(dense$logscale_design$storage, "dense")
  expect_equal(dense$logscale_design$dims, c(2L, 2L))
})

test_that("bundled subset header exposes thread-local install and clear hooks", {
  path <- pdmp_subsample_hpp_path()
  code <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(code, "thread_local")
  expect_match(code, "pdmp_set_subsample_indices", fixed = TRUE)
  expect_match(code, "pdmp_clear_subsample_indices", fixed = TRUE)
  expect_match(code, "pdmp_get_subsample_index", fixed = TRUE)
})
