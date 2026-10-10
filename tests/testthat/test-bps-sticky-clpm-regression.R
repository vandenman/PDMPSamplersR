test_that("BPS preserves frozen CLPM coordinates at the R/Julia bridge", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  # The failed production rows all used the 62-dimensional P = 4, T = 4
  # CLPM state.  Coordinates 20, 24, and 31 are representative off-diagonal
  # sticky coordinates from those rows.
  model <- tempfile(fileext = ".stan")
  on.exit(unlink(model), add = TRUE)
  writeLines(c(
    "data { int<lower=0, upper=1> prior_only; }",
    "parameters { vector[62] x; }",
    "model { x ~ normal(0, 1); }"
  ), model)
  data <- list(prior_only = 0L)
  # BridgeStan exposes vector coordinates as x.20, x.24, and x.31.
  sticky <- paste0("x.", c(20L, 24L, 31L))
  common <- list(
    path_to_stanmodel = model,
    standata = data,
    flow = "BouncyParticle",
    algorithm = "GridThinningStrategy",
    T = 2,
    grid_n = 6L,
    grid_t_max = 0.05,
    x0 = rep(0.01, 62L),
    theta0 = rep(-1, 62L),
    sticky = TRUE,
    can_stick = sticky,
    model_prior = betabernoulli(1, 2),
    parameter_prior = rep(dnorm(0), 62L),
    show_progress = FALSE,
    materialize = FALSE
  )

  independent <- do.call(pdmp_sample_from_stanmodel,
    c(common, list(seed = 20260702L)))
  expect_s3_class(independent, "pdmp_result")
  expect_equal(independent$d, 62L)
  expect_gt(independent$stats$sticky_events[[1L]], 0)

  dependent_common <- common
  dependent_common$parameter_prior <- NULL
  dependent <- do.call(pdmp_sample_from_stanmodel, c(dependent_common, list(
    slab_prior = dense_gaussian_slab(
      rep(0, length(sticky)), diag(length(sticky)), coef = sticky),
    seed = 20260703L
  )))
  expect_s3_class(dependent, "pdmp_result")
  expect_equal(dependent$d, 62L)
  expect_gt(dependent$stats$sticky_events[[1L]], 0)
})
