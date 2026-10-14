test_that("dependent slab constructors create validated specs", {
  dense <- dense_gaussian_slab(c(0, 1), diag(2), coef = c("b.x1", "b.x2"))
  expect_s3_class(dense, "dependent_slab_prior")
  expect_equal(dense$type, "dense_gaussian")
  expect_equal(dense$coef, c("b.x1", "b.x2"))

  exch <- exchangeable_gaussian_slab(u = 2, v = 0.1)
  expect_equal(exch$type, "exchangeable_gaussian")

  logscale <- loglinear_gaussian_scale_slab(c(0, 1), logscale = "log_tau",
    logscale_design = matrix(1, 2, 1), coef = c("b.x1", "b.x2"))
  expect_equal(logscale$type, "loglinear_gaussian_scale")
})

test_that("dependent slab constructors reject invalid inputs", {
  expect_error(dense_gaussian_slab(0, matrix(1, 2, 2)), "dimensions|dim")
  expect_error(dense_gaussian_slab(c(0, 0), matrix(c(1, 2, 0, 1), 2)), "symmetric")
  expect_equal(exchangeable_gaussian_slab(mean = 1, u = 1, v = 0)$mean, 1)
  expect_error(exchangeable_gaussian_slab(u = 1, v = 0, coef = list(a = 1)), "coef")
  expect_error(loglinear_gaussian_scale_slab(0, list(a = 1), matrix(1, 1, 1)), "logscale")
})

test_that("validate_pdmp_params separates legacy and dependent sticky modes", {
  d <- 3
  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      d, "ZigZag", "GridThinningStrategy", 10,
      sticky = TRUE, can_stick = c(TRUE, TRUE, FALSE),
      model_prior = bernoulli(0.5),
      parameter_prior = rep(1, d),
      slab_prior = exchangeable_gaussian_slab(u = 1, v = 0)
    ),
    "either"
  )

  params <- PDMPSamplersR:::validate_pdmp_params(
    d, "ZigZag", "GridThinningStrategy", 10,
    sticky = TRUE, can_stick = c(TRUE, TRUE, FALSE),
    model_prior = betabernoulli(1, 1),
    slab_prior = exchangeable_gaussian_slab(u = 1, v = 0)
  )
  expect_true(PDMPSamplersR:::is.slab_prior(params$slab_prior))
  expect_null(params$parameter_prior)

  adaptive <- PDMPSamplersR:::validate_pdmp_params(
    d, "AdaptiveBoomerang", "GridThinningStrategy", 10,
    sticky = TRUE, can_stick = c(TRUE, TRUE, FALSE),
    model_prior = bernoulli(0.5),
    slab_prior = exchangeable_gaussian_slab(u = 1, v = 0)
  )
  expect_identical(adaptive$flow, "AdaptiveBoomerang")
  expect_silent(PDMPSamplersR:::.validate_dependent_slab_early(
    adaptive$slab_prior, TRUE, adaptive$flow, adaptive$algorithm,
    adaptive$model_prior, adaptive$parameter_prior))
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
      sticky = TRUE, can_stick = c(TRUE, TRUE, FALSE, FALSE),
      model_prior = bernoulli(0.5),
      slab_prior = loglinear_gaussian_scale_slab(0, logscale = 1,
        logscale_design = matrix(1, 1, 1), coef = 1)
    ),
    "disjoint"
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
})

test_that("public custom-gradient dependent slab has one full-target contract", {
  expect_false("prior_grad" %in% names(formals(pdmp_sample)))
})

test_that("public custom-gradient dependent slab path can run a tiny chain", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  result <- pdmp_sample(
    function(x) x,
    d = 2,
    flow = "ZigZag",
    algorithm = "GridThinningStrategy",
    T = 1,
    sticky = TRUE,
    can_stick = c(TRUE, TRUE),
    model_prior = bernoulli(0.5),
    slab_prior = dense_gaussian_slab(c(0, 0), diag(2)),
    show_progress = FALSE,
    materialize = FALSE
  )
  expect_s3_class(result, "pdmp_result")
})

test_that("dependent slabs run through AdaptiveBoomerang warmup and sampling", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  run_chain <- function(scheme, seed, n_chains = 1L) {
    pdmp_sample(
      function(x) 2 * x,
      d = 2,
      flow = "AdaptiveBoomerang",
      algorithm = "GridThinningStrategy",
      adaptive_scheme = scheme,
      T = 30, t_warmup = 10,
      x0 = c(0.4, -0.3),
      sticky = TRUE,
      can_stick = c(TRUE, TRUE),
      model_prior = bernoulli(0.8),
      slab_prior = dense_gaussian_slab(c(0, 0), diag(2)),
      n_chains = n_chains, seed = seed,
      show_progress = FALSE, materialize = FALSE
    )
  }

  for (scheme in c("diagonal", "fullrank")) {
    result <- run_chain(scheme, 810L + match(scheme,
      c("diagonal", "fullrank")))
    expect_s3_class(result, "pdmp_result")
    expect_equal(result$d, 2L)
    numeric_stats <- unlist(Filter(is.numeric, result$stats))
    expect_true(all(is.finite(numeric_stats)))
    expect_true(all(result$stats$warmup_events > 0))
    expect_true(all(result$stats$main_events > 0))
    expect_true(all(abs(result$stats$final_lambda_ref - 0.1) > 1e-8))
    inclusion <- inclusion_probs(result)
    expect_length(inclusion, 2L)
    expect_true(all(is.finite(inclusion)))
    expect_true(all(inclusion >= 0 & inclusion <= 1))
  }

  multi1 <- run_chain("diagonal", 812L, n_chains = 2L)
  multi2 <- run_chain("diagonal", 812L, n_chains = 2L)
  # Wall-clock timings and memory use differ between identical runs.
  deterministic_names <- grep("elapsed_time|_seconds$|_allocated_bytes$",
                              names(multi1$stats), value = TRUE, invert = TRUE)
  expect_equal(multi1$stats[deterministic_names],
    multi2$stats[deterministic_names])
  expect_equal(inclusion_probs(multi1), inclusion_probs(multi2))
  expect_length(multi1$stats$final_lambda_ref, 2L)
  expect_true(all(abs(multi1$stats$final_lambda_ref - 0.1) > 1e-8))
})

test_that("Stan-backed dependent slabs validate files before Julia", {
  expect_error(
    pdmp_sample_from_stanmodel(
      "missing.stan", "missing.json",
      sticky = TRUE,
      algorithm = "GridThinningStrategy",
      model_prior = bernoulli(0.5),
      slab_prior = dense_gaussian_slab(0, matrix(1, 1, 1), coef = 1)
    ),
    "not found"
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

  logscale <- loglinear_gaussian_scale_slab(0, logscale = "b.Intercept",
    logscale_design = matrix(1, 2, 1), coef = c("b.x1", "b.x2"))
  JuliaCall::julia_assign("r_logscale_slab", logscale)
  logscale_provider_type <- JuliaCall::julia_eval("string(typeof(build_slab_provider(r_logscale_slab, r_unc_names, r_can_stick, 3)))")
  expect_match(logscale_provider_type, "LogLinearGaussianScaleSlab")
  logscale_clock_type <- JuliaCall::julia_eval("
    string(typeof(default_aggregate_unstick_clock(
      build_slab_provider(r_logscale_slab, r_unc_names, r_can_stick, 3),
      build_model_prior_odds(r_model_prior, [2, 3], 3)
    )))
  ")
  expect_match(logscale_clock_type, "ExponentialSumAggregateClock")

  if (requireNamespace("Matrix", quietly = TRUE)) {
    sparse_design <- Matrix::sparseMatrix(
      i = c(1L, 2L, 2L), j = c(1L, 1L, 2L),
      x = c(1, 0.5, 1), dims = c(2L, 2L)
    )
    structured <- loglinear_gaussian_scale_slab(
      0, logscale = c("b.Intercept", "log_lambda"),
      logscale_design = sparse_design, coef = c("b.x1", "b.x2")
    )
    JuliaCall::julia_assign("r_structured_sparse_slab", structured)
    expect_true(JuliaCall::julia_eval("
      begin
        provider = build_slab_provider(
          r_structured_sparse_slab,
          [\"b.x1\", \"b.x2\", \"b.Intercept\", \"log_lambda\"],
          BitVector([true, true, false, false]), 4
        )
        !hasfield(typeof(provider), :logscale_design) &&
          length(provider.nzval) == 3 && provider.rowptr == [1, 2, 4]
      end
    "))
  }

  alg_type <- JuliaCall::julia_eval("
    string(typeof(wrap_dependent_sticky(
      GridThinningStrategy(), true, r_model_prior, r_dense_slab,
      r_can_stick, build_flow(\"AdaptiveBoomerang\",
        Matrix{Float64}(I, 3, 3), zeros(3);
        adaptive_scheme=\"diagonal\"), r_unc_names
    )))
  ")
  expect_match(alg_type, "AggregateSticky")

  loglinear <- loglinear_gaussian_scale_slab(
    0, logscale = "b.Intercept", logscale_design = matrix(1, 2, 1),
    coef = c("b.x1", "b.x2"))
  JuliaCall::julia_assign("r_loglinear_slab", loglinear)
  adaptive_bridge <- JuliaCall::julia_eval("
    begin
      flow = build_flow(\"AdaptiveBoomerang\", Matrix{Float64}(I, 3, 3),
        zeros(3); adaptive_scheme=\"diagonal\")
      alg = wrap_dependent_sticky(GridThinningStrategy(), true,
        r_model_prior, r_loglinear_slab, r_can_stick, flow, r_unc_names)
      flow_copy = PDMPSamplers._copy_flow(flow)
      alg_copy = PDMPSamplers._copy_algorithm(alg)
      (string(typeof(flow)), string(typeof(alg.clock)), flow !== flow_copy,
       alg.clock !== alg_copy.clock &&
         alg.clock.workspace !== alg_copy.clock.workspace)
    end
  ")
  expect_match(adaptive_bridge[[1]], "MutableBoomerang")
  expect_match(adaptive_bridge[[2]], "HarmonicLogLinearAggregateClock")
  expect_true(adaptive_bridge[[3]])
  expect_true(adaptive_bridge[[4]])
})

