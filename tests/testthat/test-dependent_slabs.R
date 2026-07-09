test_that("dependent slab constructors create validated specs", {
  indep <- independent_slab_density(c(0.5, 1), coef = c("b.x1", "b.x2"))
  expect_s3_class(indep, "dependent_slab_prior")
  expect_equal(indep$type, "independent_slab_density")
  expect_equal(indep$coef, c("b.x1", "b.x2"))

  dense <- dense_gaussian_slab(c(0, 1), diag(2), coef = 2:3)
  expect_s3_class(dense, "dependent_slab_prior")
  expect_equal(dense$type, "dense_gaussian")

  exch <- exchangeable_gaussian_slab(u = 2, v = 0.1)
  expect_equal(exch$type, "exchangeable_gaussian")
  expect_true(exch$zero_mean)
})

test_that("dependent slab constructors reject invalid inputs", {
  expect_error(independent_slab_density(0), "positive")
  expect_error(dense_gaussian_slab(0, matrix(1, 2, 2)), "dimensions|dim")
  expect_error(dense_gaussian_slab(c(0, 0), matrix(c(1, 2, 0, 1), 2)), "symmetric")
  expect_error(exchangeable_gaussian_slab(mean = 1, u = 1, v = 0), "mean")
  expect_error(exchangeable_gaussian_slab(u = 1, v = 0, coef = list(a = 1)), "coef")
})

test_that("validate_pdmp_params separates legacy and dependent sticky modes", {
  d <- 3
  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      d, "ZigZag", "GridThinningStrategy", 10,
      sticky = TRUE, can_stick = c(TRUE, TRUE, FALSE),
      model_prior = bernoulli(0.5),
      parameter_prior = rep(1, d),
      slab_prior = independent_slab_density(1)
    ),
    "either"
  )

  params <- PDMPSamplersR:::validate_pdmp_params(
    d, "ZigZag", "GridThinningStrategy", 10,
    sticky = TRUE, can_stick = c(TRUE, TRUE, FALSE),
    model_prior = exchangeable_model_size_prior(c(1, 1, 1)),
    slab_prior = exchangeable_gaussian_slab(u = 1, v = 0)
  )
  expect_true(PDMPSamplersR:::is.slab_prior(params$slab_prior))
  expect_null(params$parameter_prior)

  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      d, "Boomerang", "GridThinningStrategy", 10,
      sticky = TRUE, can_stick = c(TRUE, TRUE, FALSE),
      model_prior = bernoulli(0.5),
      slab_prior = independent_slab_density(1)
    ),
    "ZigZag and BouncyParticle"
  )
})

test_that("slab_prior requires sticky and reconciles coef with can_stick", {
  d <- 4
  slab <- dense_gaussian_slab(c(0, 0), diag(2), coef = 2:3)

  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      d, "ZigZag", "GridThinningStrategy", 10,
      sticky = FALSE, model_prior = bernoulli(0.5), slab_prior = slab
    ),
    "sticky"
  )

  params <- PDMPSamplersR:::validate_pdmp_params(
    d, "ZigZag", "GridThinningStrategy", 10,
    sticky = TRUE, can_stick = NULL,
    model_prior = bernoulli(c(0.2, 0.8)),
    slab_prior = slab
  )
  expect_equal(params$can_stick, c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(params$model_prior$prob, c(0.2, 0.8))

  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      d, "ZigZag", "GridThinningStrategy", 10,
      sticky = TRUE, can_stick = c(FALSE, TRUE, FALSE, FALSE),
      model_prior = bernoulli(c(0.2, 0.8)),
      slab_prior = slab
    ),
    "subset"
  )

  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      d, "ZigZag", "GridThinningStrategy", 10,
      sticky = TRUE, can_stick = NULL,
      model_prior = bernoulli(0.5),
      slab_prior = dense_gaussian_slab(0, matrix(1, 1, 1), coef = 5)
    ),
    "1:d"
  )

  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      d, "ZigZag", "GridThinningStrategy", 10,
      sticky = TRUE, can_stick = NULL,
      model_prior = bernoulli(0.5),
      slab_prior = dense_gaussian_slab(0, matrix(1, 1, 1), coef = "b.x")
    ),
    "unconstrained parameter names"
  )

  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      d, "ZigZag", "GridThinningStrategy", 10,
      sticky = TRUE, can_stick = NULL,
      model_prior = exchangeable_model_size_prior(c(1, 1)),
      slab_prior = slab
    ),
    "beta dimension"
  )
})

test_that("public custom-gradient sampling gates slab_prior until target correction exists", {
  expect_error(
    pdmp_sample(
      function(x) x,
      d = 2,
      sticky = TRUE,
      model_prior = bernoulli(0.5),
      slab_prior = dense_gaussian_slab(c(0, 0), diag(2), coef = 1:2)
    ),
    "not yet supported"
  )
})

test_that("Stan-backed dependent slabs are gated until target correction exists", {
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", "missing.json",
      sticky = TRUE,
      model_prior = bernoulli(0.5),
      slab_prior = dense_gaussian_slab(0, matrix(1, 1, 1), coef = 1)
    ),
    "not yet supported"
  )
})

test_that("Julia bridge builds dependent slab concrete types", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  dense <- dense_gaussian_slab(c(0, 0), diag(2), coef = c("b.x1", "b.x2"))
  prior <- bernoulli(c(0.2, 0.8, 0.5))
  can_stick <- c(FALSE, TRUE, TRUE)
  unc_names <- c("b.Intercept", "b.x1", "b.x2")

  JuliaCall::julia_assign("r_dense_slab", dense)
  JuliaCall::julia_assign("r_model_prior", prior)
  JuliaCall::julia_assign("r_can_stick", can_stick)
  JuliaCall::julia_assign("r_unc_names", unc_names)
  provider_type <- JuliaCall::julia_eval("string(typeof(build_slab_provider(r_dense_slab, r_unc_names, r_can_stick, 3)))")
  expect_match(provider_type, "DenseGaussianSlab")

  alg_type <- JuliaCall::julia_eval("
    string(typeof(wrap_dependent_sticky(
      GridThinningStrategy(), true, r_model_prior, r_dense_slab,
      r_can_stick, \"ZigZag\", r_unc_names
    )))
  ")
  expect_match(alg_type, "AggregateSticky")
})
