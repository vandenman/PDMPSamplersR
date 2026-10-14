test_that("resolved Julia flows preserve requested means and covariances", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()

  d <- 5L
  mean <- seq(-0.4, 0.4, length.out = d)
  diagonal <- diag(c(0.25, 0.5, 1, 2, 4))
  identity <- diag(d)
  dense <- crossprod(matrix(c(
    1, .1, 0, 0, 0,
    0, 1, .2, 0, 0,
    0, 0, 1, .1, 0,
    0, 0, 0, 1, .2,
    0, 0, 0, 0, 1), d, d))

  cases <- list(
    list("PreconditionedZigZag", identity),
    list("PreconditionedZigZag", diagonal),
    list("PreconditionedBPS", diagonal),
    list("Boomerang", diagonal),
    list("DensePreconditionedZigZag", dense),
    list("DensePreconditionedBPS", dense)
  )
  for (case in cases) {
    answer <- PDMPSamplersR:::.pdmpsamplers_julia_call(
      "r_resolve_flow_spec", case[[1L]], mean, case[[2L]],
      adaptive_scheme = "diagonal", refresh_rate = 0.3)
    if (is.environment(answer)) answer <- as.list(answer)
    expect_equal(as.numeric(answer$mean), mean, tolerance = 1e-12)
    expect_equal(as.matrix(answer$covariance), case[[2L]], tolerance = 1e-12,
                 ignore_attr = TRUE)
  }
})

test_that("adaptive flow honors its supplied initial specification", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()
  d <- 4L
  answer <- PDMPSamplersR:::.pdmpsamplers_julia_call(
    "r_resolve_flow_spec", "AdaptiveBoomerang", rep(0.2, d),
    diag(seq_len(d)), adaptive_scheme = "diagonal", refresh_rate = 0.1)
  if (is.environment(answer)) answer <- as.list(answer)
  expect_equal(as.numeric(answer$mean), rep(0.2, d), tolerance = 1e-12)
  expect_equal(as.matrix(answer$covariance), diag(seq_len(d)),
               tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("frozen multi-chain verification checks every resolved chain", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()
  d <- 3L
  requested_mean <- c(-0.2, 0.1, 0.4)
  requested_cov <- diag(c(0.5, 1, 2))
  resolve <- function(covariance) {
    value <- PDMPSamplersR:::.pdmpsamplers_julia_call(
      "r_resolve_flow_spec", "PreconditionedBPS", requested_mean,
      covariance, adaptive_scheme = "diagonal", refresh_rate = 0.2)
    if (is.environment(value)) value <- as.list(value)
    value
  }
  agreement <- PDMPSamplersR:::.flow_spec_agreement(
    list(resolve(requested_cov), resolve(diag(c(0.5, 1, 3)))),
    requested_mean, requested_cov)
  expect_identical(agreement$chain, c(TRUE, FALSE))
})

test_that("adaptive sampling treats supplied flow values as initialization", {
  skip_on_cran()
  skip_if_no_pdmp_julia_backend()
  model_dir <- system.file("stan", "models", package = "PDMPSamplersR")
  model_path <- file.path(model_dir, "mvnormal.stan")
  skip_if_not(file.exists(model_path), "mvnormal.stan not found")
  d <- 2L
  initial_mean <- c(0.7, -0.4)
  initial_covariance <- diag(c(0.4, 2.5))
  result <- pdmp_sample_from_stanmodel(model_path,
    list(N = d, mu = c(-0.2, 0.3), sigma = diag(d)),
    flow = "AdaptiveBoomerang", algorithm = "GridThinningStrategy",
    T = 0.25, t_warmup = 0.2, warmup_adaptation_interval = 0,
    flow_mean = initial_mean, flow_cov = initial_covariance,
    flow_specification = "initial", show_progress = FALSE,
    seed = 991L, materialize = FALSE)
  provenance <- attr(result, "flow_provenance")
  expect_identical(provenance$specification_role, "initial")
  expect_true(provenance$preflight_agreement)
  expect_true(provenance$verification_passed)
  expect_identical(provenance$n_verified_chains, 1L)
  expect_length(provenance$resolved_mean_hash, 1L)
  expect_length(provenance$resolved_covariance_hash, 1L)
})
