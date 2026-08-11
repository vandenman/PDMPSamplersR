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

  logscale <- independent_logscale_gaussian_slab(log_base_scales = c(0, 1), logscale = "log_tau", coef = c("b.x1", "b.x2"))
  expect_equal(logscale$type, "independent_logscale_gaussian")

  global_exch <- global_logscale_exchangeable_gaussian_slab(logscale = "log_tau", u = 1, v = 0.2)
  expect_equal(global_exch$type, "global_logscale_exchangeable_gaussian")
})

test_that("dependent slab constructors reject invalid inputs", {
  expect_error(independent_slab_density(0), "positive")
  expect_error(dense_gaussian_slab(0, matrix(1, 2, 2)), "dimensions|dim")
  expect_error(dense_gaussian_slab(c(0, 0), matrix(c(1, 2, 0, 1), 2)), "symmetric")
  expect_error(exchangeable_gaussian_slab(mean = 1, u = 1, v = 0), "mean")
  expect_error(exchangeable_gaussian_slab(u = 1, v = 0, coef = list(a = 1)), "coef")
  expect_error(independent_logscale_gaussian_slab(0, list(a = 1)), "logscale")
  expect_error(global_logscale_exchangeable_gaussian_slab(1:2, u = 1, v = 0), "length 1")
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

  adaptive <- PDMPSamplersR:::validate_pdmp_params(
    d, "AdaptiveBoomerang", "GridThinningStrategy", 10,
    sticky = TRUE, can_stick = c(TRUE, TRUE, FALSE),
    model_prior = bernoulli(0.5),
    slab_prior = independent_slab_density(1)
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
      slab_prior = independent_logscale_gaussian_slab(0, logscale = 1, coef = 1)
    ),
    "disjoint"
  )

  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      d, "ZigZag", "GridThinningStrategy", 10,
      sticky = TRUE, can_stick = c(TRUE, FALSE, TRUE, FALSE),
      model_prior = bernoulli(0.5),
      slab_prior = global_logscale_exchangeable_gaussian_slab(logscale = 3, u = 1, v = 0, coef = 1)
    ),
    "non-stickable"
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

test_that("callback Gaussian slabs require active negative gradient for sampling", {
  callback_slab <- gaussian_scale_mixture_slab(function(x) {
    list(mean = rep(0, 2), cov = diag(2))
  })

  expect_error(
    PDMPSamplersR:::validate_pdmp_params(
      2, "ZigZag", "GridThinningStrategy", 10,
      sticky = TRUE,
      can_stick = c(TRUE, TRUE),
      model_prior = bernoulli(0.5),
      slab_prior = callback_slab
    ),
    "active_prior_neggrad"
  )

  expect_error(
    pdmp_sample(
      function(x) x,
      d = 2,
      flow = "ZigZag",
      algorithm = "GridThinningStrategy",
      sticky = TRUE,
      can_stick = c(TRUE, TRUE),
      model_prior = bernoulli(0.5),
      slab_prior = callback_slab
    ),
    "active_prior_neggrad"
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
  deterministic_names <- setdiff(names(multi1$stats),
    grep("elapsed_time", names(multi1$stats), value = TRUE))
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

  indep <- independent_slab_density(1, coef = c("b.x1", "b.x2"))
  JuliaCall::julia_assign("r_indep_slab", indep)
  indep_provider_type <- JuliaCall::julia_eval("string(typeof(build_slab_provider(r_indep_slab, r_unc_names, r_can_stick, 3)))")
  expect_match(indep_provider_type, "IndependentZeroMeanGaussianSlab")

  logscale <- independent_logscale_gaussian_slab(0, logscale = "b.Intercept", coef = c("b.x1", "b.x2"))
  JuliaCall::julia_assign("r_logscale_slab", logscale)
  logscale_provider_type <- JuliaCall::julia_eval("string(typeof(build_slab_provider(r_logscale_slab, r_unc_names, r_can_stick, 3)))")
  expect_match(logscale_provider_type, "IndependentZeroMeanLogscaleGaussianSlab")
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

  global_exch <- global_logscale_exchangeable_gaussian_slab(logscale = "b.Intercept", u = 1, v = 0.1, coef = c("b.x1", "b.x2"))
  JuliaCall::julia_assign("r_global_exch_slab", global_exch)
  global_provider_type <- JuliaCall::julia_eval("string(typeof(build_slab_provider(r_global_exch_slab, r_unc_names, r_can_stick, 3)))")
  expect_match(global_provider_type, "GlobalLogscaleExchangeableGaussianSlab")
  global_clock_type <- JuliaCall::julia_eval("
    string(typeof(default_aggregate_unstick_clock(
      build_slab_provider(r_global_exch_slab, r_unc_names, r_can_stick, 3),
      build_model_prior_odds(r_model_prior, [2, 3], 3)
    )))
  ")
  expect_match(global_clock_type, "ChebyshevResidualAggregateClock")

  alg_type <- JuliaCall::julia_eval("
    string(typeof(wrap_dependent_sticky(
      GridThinningStrategy(), true, r_model_prior, r_dense_slab,
      r_can_stick, build_flow(\"AdaptiveBoomerang\",
        Matrix{Float64}(I, 3, 3), zeros(3);
        adaptive_scheme=\"diagonal\"), r_unc_names
    )))
  ")
  expect_match(alg_type, "AggregateSticky")

  adaptive_bridge <- JuliaCall::julia_eval("
    begin
      flow = build_flow(\"AdaptiveBoomerang\", Matrix{Float64}(I, 3, 3),
        zeros(3); adaptive_scheme=\"diagonal\")
      alg = wrap_dependent_sticky(GridThinningStrategy(), true,
        r_model_prior, r_global_exch_slab, r_can_stick, flow, r_unc_names)
      flow_copy = PDMPSamplers._copy_flow(flow)
      alg_copy = PDMPSamplers._copy_algorithm(alg)
      (string(typeof(flow)), string(typeof(alg.clock)),
       alg.clock.allow_slow_fallback, flow !== flow_copy,
       alg.clock !== alg_copy.clock &&
         alg.clock.workspace !== alg_copy.clock.workspace)
    end
  ")
  expect_match(adaptive_bridge[[1]], "MutableBoomerang")
  expect_match(adaptive_bridge[[2]], "FourierResidualAggregateClock")
  expect_true(adaptive_bridge[[3]])
  expect_true(adaptive_bridge[[4]])
  expect_true(adaptive_bridge[[5]])
})

test_that("Julia bridge executes R callback slab provider contracts", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  callback_slab <- gaussian_scale_mixture_slab(
    mean_cov = function(x) {
      list(
        mean = c(0.25 + x[[2]], -0.5),
        cov = matrix(c(2.0, 0.3, 0.3, 1.5), 2, 2)
      )
    },
    active_prior_neggrad = function(x, active) {
      out <- numeric(length(x))
      if (active[[1]]) out[[1]] <- x[[1]] - 0.25
      if (active[[2]]) out[[3]] <- 2 * x[[3]]
      out
    },
    coef = c("b.one", "b.three")
  )
  arbitrary_slab <- arbitrary_slab_boundary(
    log_q_zero = function(x, active, j) {
      -0.5 * j + sum(active) + x[[2]]
    },
    active_prior_neggrad = function(x, active) {
      out <- numeric(length(x))
      if (active[[1]]) out[[1]] <- 3 * x[[1]]
      if (active[[2]]) out[[3]] <- 4 * x[[3]]
      out
    },
    coef = c("b.one", "b.three")
  )

  JuliaCall::julia_assign("r_callback_slab", callback_slab)
  JuliaCall::julia_assign("r_arbitrary_slab", arbitrary_slab)
  JuliaCall::julia_assign("r_callback_can_stick", c(TRUE, FALSE, TRUE))
  JuliaCall::julia_assign("r_callback_names", c("b.one", "b.two", "b.three"))
  JuliaCall::julia_assign("r_callback_missing_grad", gaussian_scale_mixture_slab(
    mean_cov = function(x) list(mean = c(0, 0), cov = diag(2)),
    coef = c("b.one", "b.three")
  ))

  expect_true(JuliaCall::julia_eval("
    begin
      provider = build_slab_provider(r_callback_slab, r_callback_names, r_callback_can_stick, 3)
      mean, cov = PDMPSamplers.gaussian_slab(provider, [1.0, 2.0, -1.0])
      mean ≈ [2.25, -0.5] && cov ≈ [2.0 0.3; 0.3 1.5]
    end
  "))
  expect_true(JuliaCall::julia_eval("
    begin
      provider = build_slab_provider(r_callback_slab, r_callback_names, r_callback_can_stick, 3)
      out = fill(NaN, 3)
      PDMPSamplers.active_prior_neggrad!(provider, out,
        [1.0, 2.0, -1.0], BitVector([true, true]))
      out ≈ [0.75, 0.0, -2.0]
    end
  "))
  expect_true(JuliaCall::julia_eval("
    begin
      provider = build_slab_provider(r_arbitrary_slab, r_callback_names, r_callback_can_stick, 3)
      active = BitVector([true, false])
      PDMPSamplers.log_boundary_density_zero(
        provider, [1.0, 2.0, -1.0], active, 2) ≈ 2.0
    end
  "))
  expect_true(JuliaCall::julia_eval("
    begin
      provider = build_slab_provider(r_arbitrary_slab, r_callback_names, r_callback_can_stick, 3)
      out = fill(NaN, 3)
      PDMPSamplers.active_prior_neggrad!(provider, out,
        [1.0, 2.0, -1.0], BitVector([true, true]))
      out ≈ [3.0, 0.0, -4.0]
    end
  "))
  expect_error(
    JuliaCall::julia_eval("PDMPSamplersRBridge._validate_sampling_slab_prior(r_callback_missing_grad)"),
    "active_prior_neggrad"
  )
})

test_that("R Bernoulli endpoints preserve Boomerang extended-rate semantics", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  slab <- global_logscale_exchangeable_gaussian_slab(
    logscale = "log_tau", u = 1, v = 1, logscale_offset = -10,
    coef = c("b.x1", "b.x2")
  )
  JuliaCall::julia_assign("r_endpoint_slab", slab)
  JuliaCall::julia_assign(
    "r_endpoint_names", c("b.x1", "b.x2", "log_tau"))
  JuliaCall::julia_assign("r_endpoint_mask", c(TRUE, TRUE, FALSE))

  for (endpoint in c(0, 1)) {
    JuliaCall::julia_assign("r_endpoint_prior", bernoulli(endpoint))
    for (logscale_velocity in c(0, 1)) {
      JuliaCall::julia_assign("r_endpoint_velocity", logscale_velocity)
      result <- JuliaCall::julia_eval("
        begin
          provider = build_slab_provider(r_endpoint_slab,
            r_endpoint_names, r_endpoint_mask, 3)
          odds = build_model_prior_odds(r_endpoint_prior, [1, 2], 3)
          flow = Boomerang(3)
          state = StickyPDMPState(Ref(0.0),
            SkeletonPoint([-2.0, 0.0, 0.0],
              [2.0, 0.0, r_endpoint_velocity]),
            BitVector([true, false, true]))
          mask = BitVector([true, true, false])
          clock = PDMPSamplers.FourierResidualAggregateClock(provider, odds;
            allow_slow_fallback=false)
          (PDMPSamplers.rate(clock, flow, state, 0.0, mask),
           PDMPSamplers.sample_time(Random.Xoshiro(1), clock, flow,
             state, 2.0, mask),
           PDMPSamplers.sample_time(Random.Xoshiro(1), clock, flow,
             state, Inf, mask))
        end
      ")
      if (endpoint == 0) {
        expect_equal(unlist(result), c(0, Inf, Inf))
      } else {
        expect_equal(unlist(result), c(Inf, 0, 0))
      }
    }
  }

  extreme_prior <- exchangeable_model_size_prior(c(5e-324, 5e-324, 1))
  JuliaCall::julia_assign("r_extreme_size_prior", extreme_prior)
  extreme_result <- JuliaCall::julia_eval("
    begin
      provider = build_slab_provider(r_endpoint_slab,
        r_endpoint_names, r_endpoint_mask, 3)
      odds = build_model_prior_odds(r_extreme_size_prior, [1, 2], 3)
      flow = Boomerang(3)
      state = StickyPDMPState(Ref(0.0),
        SkeletonPoint([-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]),
        BitVector([true, false, true]))
      mask = BitVector([true, true, false])
      clock = PDMPSamplers.FourierResidualAggregateClock(provider, odds;
        allow_slow_fallback=false)
      peak = PDMPSamplers._boomerang_narrow_peak_segment(
        clock, flow, state, mask)
      event_time = PDMPSamplers.sample_time(Random.Xoshiro(1), clock,
        flow, state, Inf, mask)
      (string(typeof(odds)), all(isfinite, odds.log_omega),
       peak.log_total_weight, event_time)
    end
  ")
  expect_match(extreme_result[[1]], "ExchangeableModelSizePrior")
  expect_true(extreme_result[[2]])
  expect_true(is.finite(extreme_result[[3]]))
  expect_true(is.finite(extreme_result[[4]]))
  expect_gt(extreme_result[[4]], 0)
})
