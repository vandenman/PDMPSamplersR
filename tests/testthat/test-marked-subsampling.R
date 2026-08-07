test_that("an enabled marked bank requires anchor updates", {
  data <- data.frame(y = c(0L, 1L), x = c(-1, 1))
  expect_error(
    PDMPSamplersR::brm_pdmp(
      y ~ x, data, family = brms::bernoulli(), subsample_size = 1L,
      use_anchor_bank = TRUE, n_anchor_updates = 0L,
      T = 1, t_warmup = 0.2, show_progress = FALSE
    ),
    "requires a positive.*n_anchor_updates"
  )
})

test_that("marked integer controls are validated before coercion", {
  data <- data.frame(y = c(0L, 1L), x = c(-1, 1))
  for (value in list(0, -1, 0.5, NA_real_, c(1, 2))) {
    expect_error(
      brm_pdmp(y ~ x, data, family = brms::bernoulli(),
               subsample_size = value, T = 1, show_progress = FALSE),
      "subsample_size.*positive integerish scalar"
    )
  }
  for (value in list(0, -1, 1.5, NA_real_, c(1, 2))) {
    expect_error(
      brm_pdmp(y ~ x, data, family = brms::bernoulli(),
               subsample_size = 1L, bank_capacity = value,
               T = 1, show_progress = FALSE),
      "bank_capacity.*positive integer"
    )
  }
  for (value in list(-1, 1.5, NA_real_, c(0, 1))) {
    expect_error(
      brm_pdmp(y ~ x, data, family = brms::bernoulli(),
               subsample_size = 1L, n_anchor_updates = value,
               T = 1, show_progress = FALSE),
      "n_anchor_updates.*non-negative integer"
    )
  }
})

test_that("analytic marked HCV rejects uncertified variants", {
  data <- data.frame(y = c(0L, 1L), count = c(1L, 2L), x = c(-1, 1))
  expect_error(
    PDMPSamplersR::brm_pdmp(
      count ~ x, data, family = brms::brmsfamily("poisson", "log"),
      subsample_size = 1L,
      use_hcv = TRUE, T = 1, t_warmup = 0.2, show_progress = FALSE
    ),
    "requires an affine Bernoulli-logit or binomial-logit"
  )
})

test_that("R thinning c0 remains a global bound", {
  bridge <- paste(readLines(
    system.file("julia", "main_interface_function.jl", package = "PDMPSamplersR"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(bridge, "ThinningStrategy\\(GlobalBounds\\(c0, d\\)\\)")
  expect_no_match(bridge, "GlobalBounds\\(c0 / d, d\\)")
})

test_that("marked eligibility is capability based and custom U0 is exact", {
  family <- brms::bernoulli(link = "logit")
  accepted <- list(
    y ~ 1, y ~ 0 + x, y ~ x, y ~ x + z,
    y ~ x * z, y ~ factor_group, y ~ x + offset(o)
  )
  for (formula in accepted) {
    expect_true(PDMPSamplersR:::marked_subsampling_eligibility(
      formula, family)$eligible, info = deparse(formula))
  }

  modifiers <- list(
    weights = y | weights(w) ~ x,
    subset = y | subset(use) ~ x
  )
  for (name in names(modifiers)) {
    expect_true(PDMPSamplersR:::marked_subsampling_eligibility(
      modifiers[[name]], family)$eligible, info = name)
  }

  rejected <- list(
    monotonic = y ~ x + mo(z),
    measurement_error = y ~ me(x, se),
    missing_predictor = y ~ mi(x),
    spline = y ~ s(x),
    group_effect = y ~ x + (1 | g),
    gaussian_process = y ~ gp(x),
    autocorrelation = y ~ x + ar()
  )
  for (name in names(rejected)) {
    expect_false(PDMPSamplersR:::marked_subsampling_eligibility(
      rejected[[name]], family)$eligible, info = name)
  }
  expect_true(PDMPSamplersR:::marked_subsampling_eligibility(
    y | trials(n) ~ x, brms::brmsfamily("binomial", "logit"))$eligible)
  expect_true(PDMPSamplersR:::marked_subsampling_eligibility(
    y | rate(exposure) ~ x, brms::brmsfamily("poisson", "log"))$eligible)
  expect_true(PDMPSamplersR:::marked_subsampling_eligibility(
    y ~ x, brms::brmsfamily("gaussian", "identity"))$eligible)

  nonlinear <- brms::bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE)
  expect_false(PDMPSamplersR:::marked_subsampling_eligibility(
    nonlinear, family)$eligible)
  expect_false(PDMPSamplersR:::marked_subsampling_eligibility(
    y | cens(censoring) ~ x, brms::brmsfamily("gaussian", "identity"))$eligible)
  expect_false(PDMPSamplersR:::marked_subsampling_eligibility(
    y | trunc(lb = 0) ~ x, brms::brmsfamily("gaussian", "identity"))$eligible)
  unsupported_family <- brms::custom_family(
    "opaque_family", dpars = "mu", links = "identity",
    type = "real", vars = "vreal1[n]"
  )
  expect_false(PDMPSamplersR:::marked_subsampling_eligibility(
    y ~ x, unsupported_family)$eligible)

  observation_code <- brms::stanvar(
    scode = "target += normal_lpdf(Y | rep_vector(0, N), 1);",
    block = "model"
  )
  expect_true(PDMPSamplersR:::marked_subsampling_eligibility(
    y ~ x, family, observation_code)$eligible)
  opaque_global_code <- brms::stanvar(
    scode = "target += normal_lpdf(auxiliary_parameter | 0, 1);",
    block = "model"
  )
  expect_true(PDMPSamplersR:::marked_subsampling_eligibility(
    y ~ x, family, opaque_global_code)$eligible)

  custom_observations <- brms::stanvar(1:4, name = "custom_y") +
    brms::stanvar(
      scode = "target += normal_lpdf(custom_y | rep_vector(0, N), 1);",
      block = "model"
    )
  expect_true(PDMPSamplersR:::marked_subsampling_eligibility(
    y ~ x, family, custom_observations, list(N = 4L))$eligible)

  conditional <- brms::stanvar(
    scode = "if (!prior_only) target += normal_lpdf(b[1] | 0, 0.5);",
    block = "model"
  )
  result <- PDMPSamplersR:::marked_subsampling_eligibility(
    y ~ x, family, conditional)
  expect_false(result$eligible)
  expect_match(result$reason, "prior_only")

  predictor_change <- brms::stanvar(
    scode = "mu += rep_vector(0.1, N);", block = "model"
  )
  expect_false(PDMPSamplersR:::marked_subsampling_eligibility(
    y ~ x, family, predictor_change)$eligible)
})

test_that("response modifiers produce exact nonnegative multipliers", {
  expect_equal(PDMPSamplersR:::marked_observation_multipliers(
    list(N = 3L, weights = c(0.5, 2, 1)), "bernoulli"), c(0.5, 2, 1))
  expect_equal(PDMPSamplersR:::marked_observation_multipliers(
    list(N = 3L, weights = c(0.5, 2, 1), trials = c(2, 3, 4)),
    "binomial"), c(1, 6, 4))
  expect_error(PDMPSamplersR:::marked_observation_multipliers(
    list(N = 2L, weights = c(1, -1)), "bernoulli"), "nonnegative")
})

test_that("centered brms design maps to unconstrained coordinates", {
  X <- cbind(Intercept = 1, x = c(-2, 0, 5), interaction = c(3, 1, -1))
  attr(X, "assign") <- c(0L, 1L, 2L)
  sdata <- list(N = 3L, X = X, means_X = c(1, 1))
  names <- c("b.1", "unused_prior_parameter", "Intercept", "b.2")
  predictor <- PDMPSamplersR:::build_marked_compact_affine_design(
    sdata, names, TRUE, "mu"
  )

  expect_equal(predictor$indices, c(1L, 3L, 4L))
  expect_equal(predictor$design[, 1], X[, 2] - 1)
  expect_equal(predictor$design[, 2], rep(1, 3))
  expect_equal(predictor$design[, 3], X[, 3] - 1)

  X_no_intercept <- cbind(x = c(-2, 0, 5), interaction = c(3, 1, -1))
  attr(X_no_intercept, "assign") <- c(1L, 2L)
  predictor_no_intercept <- PDMPSamplersR:::build_marked_compact_affine_design(
    list(N = 3L, X = X_no_intercept), c("b.1", "b.2"), FALSE, "mu"
  )
  expect_equal(predictor_no_intercept$indices, c(1L, 2L))
  expect_equal(dim(predictor_no_intercept$design), dim(X_no_intercept))
  expect_equal(as.numeric(predictor_no_intercept$design), as.numeric(X_no_intercept))
})

test_that("opaque deterministic data are not reduced to one observation", {
  sdata <- list(
    N = 4L, Y = c(0L, 1L, 1L, 0L), X = matrix(1:8, 4, 2),
    custom_global_data = c(2, 3, 5, 7), prior_only = 0L
  )
  prior <- PDMPSamplersR:::make_opaque_deterministic_standata(sdata)
  expect_equal(prior$N, sdata$N)
  expect_identical(prior$Y, sdata$Y)
  expect_identical(prior$custom_global_data, sdata$custom_global_data)
  expect_equal(prior$prior_only, 1L)
})

test_that("Julia-owned N over m scaling closes the full gradient", {
  Z <- cbind(1, c(-1.5, -0.5, 0.5, 2))
  y <- c(0, 1, 0, 1)
  offset <- c(0.2, -0.1, 0.3, 0)
  anchor <- c(-0.2, 0.1)
  theta <- c(0.4, -0.3)
  prior <- function(x) c(x[1] / 4, x[2] / 9)
  likelihood_grad <- function(x, idx) {
    eta <- drop(Z[idx, , drop = FALSE] %*% x) + offset[idx]
    colSums(Z[idx, , drop = FALSE] * (plogis(eta) - y[idx]))
  }
  full <- prior(theta) + likelihood_grad(theta, seq_len(nrow(Z)))
  full_anchor <- prior(anchor) + likelihood_grad(anchor, seq_len(nrow(Z)))
  deterministic <- full_anchor + prior(theta) - prior(anchor)

  subsets <- combn(seq_len(nrow(Z)), 2L, simplify = FALSE)
  estimates <- vapply(subsets, function(S) {
    residual <- likelihood_grad(theta, S) - likelihood_grad(anchor, S)
    deterministic + nrow(Z) / length(S) * residual
  }, numeric(2))
  expect_equal(rowMeans(estimates), full, tolerance = 1e-12)
})

test_that("generated brms affine predictors close and satisfy the envelope", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("PDMPSAMPLERSR_BRMS_CLOSURE_TESTS"), "true"),
    "Generated BridgeStan closure tests are disabled"
  )
  skip_if_no_brms_setup()

  data <- data.frame(
    y = c(0L, 1L, 0L, 1L, 1L, 0L),
    x = c(-1.2, -0.5, 0.1, 0.4, 1.1, 1.7),
    z = c(-0.7, 0.2, 1.0, -1.1, 0.5, 1.4),
    factor_group = factor(c("a", "b", "c", "a", "c", "b")),
    o = c(0.2, -0.1, 0.3, 0, -0.2, 0.1)
  )
  formulas <- list(
    intercept = y ~ 1,
    no_intercept = y ~ 0 + x,
    fixed = y ~ x,
    multiple = y ~ x + z,
    interaction = y ~ x * z,
    factor = y ~ factor_group,
    offset = y ~ x + offset(o)
  )
  family <- brms::bernoulli(link = "logit")

  for (name in names(formulas)) {
    formula <- formulas[[name]]
    eligibility <- PDMPSamplersR:::marked_subsampling_eligibility(
      formula, family
    )
    expect_true(eligibility$eligible, info = name)
    scode <- brms::stancode(formula, data = data, family = family)
    sdata <- brms::standata(formula, data = data, family = family)
    prior_data <- PDMPSamplersR:::make_opaque_deterministic_standata(sdata)
    stan_file <- PDMPSamplersR:::cached_stan_model(scode)
    full_file <- tempfile(fileext = ".json")
    prior_file <- tempfile(fileext = ".json")
    PDMPSamplersR:::write_stan_json(sdata, full_file)
    PDMPSamplersR:::write_stan_json(prior_data, prior_file)
    unc_names <- PDMPSamplersR:::.pdmpsamplers_julia_call(
      "r_get_param_unc_names", normalizePath(stan_file), normalizePath(full_file)
    )
    geometry <- PDMPSamplersR:::build_marked_predictor_geometry(
      sdata, unc_names, eligibility
    )
    anchor <- seq(-0.15, 0.05, length.out = length(unc_names))
    theta <- seq(0.1, 0.3, length.out = length(unc_names))
    diagnostics <- PDMPSamplersR:::.pdmpsamplers_julia_call(
      "r_marked_family_closure_diagnostics",
      normalizePath(stan_file), normalizePath(full_file), normalizePath(prior_file),
      eligibility$family, geometry$designs, geometry$design_indices,
      geometry$dimension, geometry$offsets, geometry$response, geometry$se,
      rep(1, nrow(data)), anchor, theta, FALSE
    )
    expect_lt(diagnostics$closure_error, 1e-8)
    expect_lte(diagnostics$max_bound_ratio, 1 + 1e-10)
    expect_equal(diagnostics$analytic_residual_calls, nrow(data))
  }

  conditional <- brms::stanvar(
    scode = "if (!prior_only) target += normal_lpdf(b[1] | 0, 0.5);",
    block = "model"
  )
  formula <- y ~ x
  eligibility <- PDMPSamplersR:::marked_subsampling_eligibility(
    formula, family, conditional)
  expect_false(eligibility$eligible)
  scode <- brms::stancode(formula, data = data, family = family,
                          stanvars = conditional)
  sdata <- brms::standata(formula, data = data, family = family,
                          stanvars = conditional)
  stan_file <- PDMPSamplersR:::cached_stan_model(scode)
  full_file <- tempfile(fileext = ".json"); prior_file <- tempfile(fileext = ".json")
  PDMPSamplersR:::write_stan_json(sdata, full_file)
  PDMPSamplersR:::write_stan_json(
    PDMPSamplersR:::make_opaque_deterministic_standata(sdata), prior_file)
  unc_names <- PDMPSamplersR:::.pdmpsamplers_julia_call(
    "r_get_param_unc_names", normalizePath(stan_file), normalizePath(full_file))
  predictor <- PDMPSamplersR:::build_marked_compact_affine_design(
    sdata, unc_names, TRUE, "mu"
  )
  diagnostics <- PDMPSamplersR:::.pdmpsamplers_julia_call(
    "r_marked_family_closure_diagnostics", normalizePath(stan_file),
    normalizePath(full_file), normalizePath(prior_file), "bernoulli",
    list(predictor$design), list(predictor$indices), length(unc_names),
    matrix(0, nrow(data), 1L), matrix(sdata$Y, nrow(data), 1L), numeric(),
    rep(1, nrow(data)), rep(0, length(unc_names)),
    rep(0.1, length(unc_names)), FALSE)
  expect_gt(diagnostics$closure_error, 0.1)
})

test_that("first family batch closes against generated brms models", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("PDMPSAMPLERSR_BRMS_CLOSURE_TESTS"), "true"),
    "Generated BridgeStan closure tests are disabled"
  )
  skip_if_no_brms_setup()
  data <- data.frame(
    y = c(1, 2, 0, 1), x = c(-1, 0.2, 0.7, 1.3),
    yb = c(1L, 1L, 0L, 1L),
    z = c(0.1, -0.2, 0.3, 0.4), w = c(0.5, 2, 1, 1),
    n = c(2, 3, 2, 4), exposure = 1:4, sev = c(0.5, 0.7, 0.6, 0.8),
    os = c(-0.3, 0.1, 0.25, -0.15)
  )
  cases <- list(
    bernoulli = list(yb ~ x,
                     brms::brmsfamily("bernoulli", "logit"), data),
    binomial = list(y | trials(n) + weights(w) ~ x,
                    brms::brmsfamily("binomial", "logit"), data),
    poisson_rate = list(y | rate(exposure) ~ x,
                        brms::brmsfamily("poisson", "log"), data),
    gaussian = list(y ~ x, brms::brmsfamily("gaussian", "identity"), data),
    gaussian_distributional = list(brms::bf(y ~ x, sigma ~ z),
                                   brms::brmsfamily("gaussian", "identity"), data),
    gaussian_distributional_sigma_offset = list(
      brms::bf(y ~ x, sigma ~ z + offset(os)),
      brms::brmsfamily("gaussian", "identity"), data),
    gaussian_known_se = list(y | se(sev, sigma = FALSE) ~ x,
                             brms::brmsfamily("gaussian", "identity"), data),
    gaussian_se_unknown_scale = list(y | se(sev, sigma = TRUE) ~ x,
                                     brms::brmsfamily("gaussian", "identity"), data),
    categorical = list(
      y ~ x, brms::brmsfamily("categorical", "logit"),
      data.frame(y = factor(c("a", "b", "c", "a")), x = data$x)
    ),
    multinomial = list(
      cbind(y1, y2, y3) | trials(n) ~ x,
      brms::brmsfamily("multinomial", "logit"),
      data.frame(y1 = c(1, 0, 2, 0), y2 = c(0, 2, 0, 1),
                 y3 = c(1, 0, 0, 1), n = 2, x = data$x)
    ),
    independent_multivariate = list(
      brms::bf(y1 ~ x, family = brms::bernoulli()) +
        brms::bf(y2 ~ z, family = brms::bernoulli()) + brms::set_rescor(FALSE),
      brms::brmsfamily("gaussian", "identity"),
      data.frame(y1 = c(0, 1, 0, 1), y2 = c(1, 1, 0, 0),
                 x = data$x, z = data$z)
    )
  )
  for (name in names(cases)) {
    formula <- cases[[name]][[1L]]; family <- cases[[name]][[2L]]
    case_data <- cases[[name]][[3L]]
    sdata <- brms::standata(formula, case_data, family = family)
    eligibility <- PDMPSamplersR:::marked_subsampling_eligibility(
      formula, family, sdata = sdata)
    expect_true(eligibility$eligible, info = name)
    stan_file <- PDMPSamplersR:::cached_stan_model(
      brms::stancode(formula, case_data, family = family))
    full_file <- tempfile(fileext = ".json"); prior_file <- tempfile(fileext = ".json")
    PDMPSamplersR:::write_stan_json(sdata, full_file)
    PDMPSamplersR:::write_stan_json(
      PDMPSamplersR:::make_opaque_deterministic_standata(sdata), prior_file)
    names_unc <- PDMPSamplersR:::.pdmpsamplers_julia_call(
      "r_get_param_unc_names", normalizePath(stan_file), normalizePath(full_file))
    geometry <- PDMPSamplersR:::build_marked_predictor_geometry(
      sdata, names_unc, eligibility)
    expect_type(geometry$designs, "list")
    expect_equal(length(geometry$designs), length(geometry$design_indices),
                 info = name)
    expect_true(all(vapply(seq_along(geometry$designs), function(k) {
      nrow(geometry$designs[[k]]) == sdata$N &&
        ncol(geometry$designs[[k]]) == length(geometry$design_indices[[k]])
    }, logical(1))), info = name)
    multipliers <- PDMPSamplersR:::marked_observation_multipliers(
      sdata, eligibility$family)
    anchors <- list(
      rep(0, length(names_unc)),
      seq(-0.03, 0.03, length.out = length(names_unc))
    )
    for (anchor in anchors) {
      diagnostics <- PDMPSamplersR:::.pdmpsamplers_julia_call(
        "r_marked_family_closure_diagnostics", normalizePath(stan_file),
        normalizePath(full_file), normalizePath(prior_file), eligibility$family,
        geometry$designs, geometry$design_indices, as.integer(geometry$dimension),
        geometry$offsets, geometry$response, geometry$se,
        multipliers, anchor, anchor + 0.05)
      expect_lt(diagnostics$closure_error, 1e-8)
      expect_lte(diagnostics$max_bound_ratio, 1 + 1e-10)
      if (eligibility$family %in% c("bernoulli", "binomial")) {
        hcv_diagnostics <- PDMPSamplersR:::.pdmpsamplers_julia_call(
          "r_marked_family_closure_diagnostics", normalizePath(stan_file),
          normalizePath(full_file), normalizePath(prior_file), eligibility$family,
          geometry$designs, geometry$design_indices,
          as.integer(geometry$dimension), geometry$offsets,
          geometry$response, geometry$se, multipliers, anchor, anchor + 0.05,
          TRUE
        )
        expect_lt(hcv_diagnostics$closure_error, 1e-8)
        expect_lte(hcv_diagnostics$max_bound_ratio, 1 + 1e-10)
      }
    }
  }
})

test_that("weighted Bernoulli models use marked acceleration", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("PDMPSAMPLERSR_SLOW_TESTS"), "true"),
    "Slow BridgeStan fallback tests are disabled"
  )
  skip_if_no_brms_setup()

  data <- data.frame(
    y = c(0L, 1L, 0L, 1L, 1L, 0L),
    x = c(-1, -0.5, 0, 0.5, 1, 1.5),
    w = c(1, 2, 1, 0.5, 1.5, 1)
  )
  fit <- PDMPSamplersR::brm_pdmp(
    y | weights(w) ~ x, data = data, family = brms::bernoulli(),
    flow = "BouncyParticle", T = 1, t_warmup = 0,
    subsample_size = 2L, show_progress = FALSE, seed = 914
  )
  expect_true(isTRUE(attr(fit, "marked_subsampling")))
})

test_that("all brms GridThinning dynamics use the marked provider", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("PDMPSAMPLERSR_SLOW_TESTS"), "true"),
              "Slow production flow tests are disabled")
  skip_if_no_brms_setup()
  data <- data.frame(y = rep(c(0L, 1L), 15), x = seq(-1, 1, length.out = 30))
  flows <- c(
    "BouncyParticle", "ZigZag", "PreconditionedBPS", "PreconditionedZigZag",
    "DensePreconditionedBPS", "DensePreconditionedZigZag",
    "Boomerang", "AdaptiveBoomerang"
  )
  for (flow in flows) {
    for (analytic_hcv in c(FALSE, TRUE)) {
      fit <- PDMPSamplersR::brm_pdmp(
        y ~ x, data, family = brms::bernoulli(), flow = flow,
        T = 0.2, t_warmup = 0.04, subsample_size = 5L,
        n_anchor_updates = 1L, use_anchor_bank = TRUE, bank_capacity = 2L,
        use_hcv = analytic_hcv,
        show_progress = FALSE, seed = 711
      )
      expect_true(isTRUE(attr(fit, "marked_subsampling")))
      expect_identical(
        attr(fit, "bridge_call_counts")[[1L]]$analytic_hcv,
        analytic_hcv
      )
    }
  }
  sticky_flows <- c("BouncyParticle", "ZigZag", "PreconditionedBPS",
                    "PreconditionedZigZag", "DensePreconditionedBPS",
                    "DensePreconditionedZigZag", "Boomerang", "AdaptiveBoomerang")
  sticky_data <- transform(
    data,
    x1 = x,
    x2 = sin(seq_along(x)),
    x3 = cos(seq_along(x)),
    x4 = seq(-1, 1, length.out = length(x))
  )
  sticky_formula <- y ~ x1 + x2 + x3 + x4
  sticky_initial <- c(0.05, -0.05, 0.05, -0.05, 0)
  for (flow_index in seq_along(sticky_flows)) {
    flow <- sticky_flows[[flow_index]]
    fit <- PDMPSamplersR::brm_pdmp(
      sticky_formula, sticky_data, family = brms::bernoulli(), flow = flow,
      prior = brms::prior(normal(0, 2), class = b),
      T = 30, t_warmup = 6, subsample_size = 5L,
      flow_mean = sticky_initial,
      use_hcv = TRUE,
      sticky = TRUE, model_prior = PDMPSamplersR::bernoulli(0.99),
      show_progress = FALSE, seed = 712 + flow_index
    )
    expect_true(isTRUE(attr(fit, "marked_subsampling")))
    expect_true(attr(fit, "pdmp_stats")$sticky_freezes >= 1, info = flow)
    expect_true(attr(fit, "pdmp_stats")$sticky_unfreezes >= 1, info = flow)
  }
  lowrank_fit <- PDMPSamplersR::brm_pdmp(
    sticky_formula, sticky_data, family = brms::bernoulli(),
    flow = "AdaptiveBoomerang", flow_mean = sticky_initial,
    prior = brms::prior(normal(0, 2), class = b),
    adaptive_scheme = "lowrank", T = 30, t_warmup = 6,
    grid_n = 100, subsample_size = 5L, use_hcv = TRUE,
    sticky = TRUE, model_prior = PDMPSamplersR::bernoulli(0.99),
    show_progress = FALSE, seed = 714
  )
  expect_true(isTRUE(attr(lowrank_fit, "marked_subsampling")))
  expect_true(attr(lowrank_fit, "pdmp_stats")$sticky_freezes >= 1,
              info = "AdaptiveBoomerang lowrank")
  expect_true(attr(lowrank_fit, "pdmp_stats")$sticky_unfreezes >= 1,
              info = "AdaptiveBoomerang lowrank")
})

test_that("all brms ThinningStrategy dynamics use the marked provider", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("PDMPSAMPLERSR_SLOW_TESTS"), "true"),
              "Slow production flow tests are disabled")
  skip_if_no_brms_setup()
  data <- data.frame(y = rep(c(0L, 1L), 15), x = seq(-1, 1, length.out = 30))
  flows <- c(
    "BouncyParticle", "ZigZag", "PreconditionedBPS", "PreconditionedZigZag",
    "DensePreconditionedBPS", "DensePreconditionedZigZag",
    "Boomerang", "AdaptiveBoomerang"
  )
  for (flow in flows) {
    fit <- PDMPSamplersR::brm_pdmp(
      y ~ x, data, family = brms::bernoulli(), flow = flow,
      algorithm = "ThinningStrategy", c0 = 100,
      T = 0.2, t_warmup = 0.04, subsample_size = 5L,
      n_anchor_updates = 1L, use_anchor_bank = TRUE, bank_capacity = 2L,
      use_hcv = TRUE, show_progress = FALSE, seed = 713
    )
    expect_true(isTRUE(attr(fit, "marked_subsampling")))
    counts <- attr(fit, "bridge_call_counts")[[1L]]
    expect_true(counts$analytic_hcv)
    expect_equal(counts$full_gradient_calls, 1 + counts$anchor_preparations)
  }

  poisson_data <- data.frame(y = rep(0:2, 10), x = seq(-1, 1, length.out = 30))
  fallback <- PDMPSamplersR::brm_pdmp(
    y ~ x, poisson_data, family = brms::brmsfamily("poisson", "log"),
    flow = "BouncyParticle", algorithm = "ThinningStrategy", c0 = 100,
    T = 0.1, t_warmup = 0, subsample_size = 5L,
    show_progress = FALSE, seed = 714
  )
  expect_false(isTRUE(attr(fallback, "marked_subsampling")))

  gaussian_data <- transform(data, y = seq(-1, 1, length.out = nrow(data)),
                             se = rep(0.7, nrow(data)))
  fixed_gaussian <- PDMPSamplersR::brm_pdmp(
    y | se(se, sigma = FALSE) ~ x, gaussian_data,
    family = brms::brmsfamily("gaussian", "identity"),
    flow = "BouncyParticle", algorithm = "ThinningStrategy", c0 = 100,
    T = 0.1, t_warmup = 0, subsample_size = 5L,
    show_progress = FALSE, seed = 715
  )
  expect_true(isTRUE(attr(fixed_gaussian, "marked_subsampling")))

  distributional_gaussian <- PDMPSamplersR::brm_pdmp(
    y ~ x, gaussian_data,
    family = brms::brmsfamily("gaussian", "identity"),
    flow = "BouncyParticle", algorithm = "ThinningStrategy", c0 = 100,
    T = 0.1, t_warmup = 0, subsample_size = 5L,
    show_progress = FALSE, seed = 716
  )
  expect_false(isTRUE(attr(distributional_gaussian, "marked_subsampling")))
})

test_that("marked anchor banks prepare during warmup and isolate chains", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("PDMPSAMPLERSR_SLOW_TESTS"), "true"),
              "Slow marked anchor-bank tests are disabled")
  skip_if_no_brms_setup()
  set.seed(803)
  data <- data.frame(x = rnorm(100))
  data$y <- rbinom(100, 1, plogis(0.4 + 0.8 * data$x))
  fit <- PDMPSamplersR::brm_pdmp(
    y ~ x, data, family = brms::bernoulli(), flow = "BouncyParticle",
    T = 4, t_warmup = 2, subsample_size = 10L,
    n_anchor_updates = 2L, use_anchor_bank = TRUE, bank_capacity = 3L,
    use_hcv = TRUE,
    n_chains = 2L, threaded = TRUE, show_progress = FALSE, seed = 804
  )
  counts <- attr(fit, "bridge_call_counts")
  expect_true(isTRUE(attr(fit, "marked_subsampling")))
  expect_length(counts, 2L)
  expect_true(all(vapply(counts, `[[`, numeric(1), "anchor_preparations") >= 1))
  expect_true(all(vapply(counts, `[[`, numeric(1), "anchor_activations") >= 1))
  expect_gte(sum(vapply(counts, `[[`, numeric(1), "anchor_main_activations")), 1)
  expect_true(all(vapply(counts, function(x) {
    x$full_gradient_calls == 1 + x$anchor_preparations
  }, logical(1))))
})

test_that("single marked anchor updates without enabling a bank", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("PDMPSAMPLERSR_SLOW_TESTS"), "true"),
              "Slow marked single-anchor tests are disabled")
  skip_if_no_brms_setup()
  data <- data.frame(y = rep(c(0L, 1L), 50), x = seq(-2, 2, length.out = 100))
  fit <- PDMPSamplersR::brm_pdmp(
    y ~ x, data, family = brms::bernoulli(), flow = "ZigZag",
    T = 4, t_warmup = 2, subsample_size = 10L,
    n_anchor_updates = 2L, use_anchor_bank = FALSE,
    use_hcv = TRUE,
    show_progress = FALSE, seed = 806
  )
  counts <- attr(fit, "bridge_call_counts")[[1L]]
  expect_true(isTRUE(attr(fit, "marked_subsampling")))
  expect_gte(counts$anchor_preparations, 1)
  expect_true(counts$analytic_hcv)
  expect_equal(counts$anchor_bank_entries, 1)
  expect_equal(counts$full_gradient_calls, 1 + counts$anchor_preparations)
})

test_that("sticky marked sampling supports anchor banks", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("PDMPSAMPLERSR_SLOW_TESTS"), "true"),
              "Slow sticky marked anchor-bank tests are disabled")
  skip_if_no_brms_setup()
  data <- data.frame(y = rep(c(0L, 1L), 50), x = seq(-2, 2, length.out = 100))
  fit <- PDMPSamplersR::brm_pdmp(
    y ~ x, data, family = brms::bernoulli(), flow = "BouncyParticle",
    prior = brms::prior(normal(0, 2), class = b),
    T = 4, t_warmup = 2, subsample_size = 10L,
    n_anchor_updates = 2L, use_anchor_bank = TRUE, bank_capacity = 3L,
    use_hcv = TRUE,
    sticky = TRUE, model_prior = PDMPSamplersR::bernoulli(0.5),
    show_progress = FALSE, seed = 807
  )
  counts <- attr(fit, "bridge_call_counts")[[1L]]
  expect_true(isTRUE(attr(fit, "marked_subsampling")))
  expect_gte(counts$anchor_preparations, 1)
  expect_true(counts$analytic_hcv)
  expect_equal(counts$full_gradient_calls, 1 + counts$anchor_preparations)
})

test_that("marked and full brms posterior moments agree", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("PDMPSAMPLERSR_SLOW_TESTS"), "true"),
    "Slow MCMC tests are disabled"
  )
  skip_if_no_brms_setup()

  set.seed(941)
  N <- 250L
  data <- data.frame(x = rnorm(N), offset = runif(N, -0.2, 0.2))
  data$y <- rbinom(N, 1, plogis(-0.35 + 0.8 * data$x + data$offset))
  # Deliberately use the optimized bernoulli_logit_glm code path. Offsets are
  # covered separately by the generated-model closure tests below.
  formula <- y ~ x
  anchor_fit <- glm(formula, data = data, family = binomial())
  anchor_coef <- coef(anchor_fit)
  marked_anchor <- c(
    unname(anchor_coef[["x"]]),
    unname(anchor_coef[["(Intercept)"]]) + mean(data$x) * anchor_coef[["x"]]
  )

  for (flow in c("BouncyParticle", "ZigZag")) {
    full <- PDMPSamplersR::brm_pdmp(
      formula, data = data, family = brms::bernoulli(), flow = flow,
      T = 2500, t_warmup = 500, show_progress = FALSE, seed = 83
    )
    marked <- PDMPSamplersR::brm_pdmp(
      formula, data = data, family = brms::bernoulli(), flow = flow,
      T = 2500, t_warmup = 500, subsample_size = 25L,
      flow_mean = marked_anchor,
      n_anchor_updates = 8L, use_anchor_bank = TRUE, bank_capacity = 4L,
      use_hcv = TRUE,
      show_progress = FALSE, seed = 83
    )
    expect_equal(
      brms::fixef(marked)[, "Estimate"],
      brms::fixef(full)[, "Estimate"],
      tolerance = 0.25
    )
    expect_true(isTRUE(attr(marked, "marked_subsampling")))
    counts <- attr(marked, "bridge_call_counts")[[1L]]
    expect_true(counts$analytic_hcv)
    expect_gte(counts$anchor_preparations, 2)
    expect_gte(counts$anchor_main_activations, 2)
    expect_equal(counts$full_gradient_calls, 1 + counts$anchor_preparations)
  }
})
