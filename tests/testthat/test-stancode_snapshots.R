# Snapshot tests for Stan code rewriting functions.
#
# These tests are gated by PDMPSAMPLERSR_EXHAUSTIVE_TESTS because they depend
# on brms code generation, which may change across versions.  Keeping them off
# CRAN avoids breakage when brms updates its Stan output.
#
# Run locally:
#   NOT_CRAN=true PDMPSAMPLERSR_EXHAUSTIVE_TESTS=true Rscript -e 'testthat::test_file("tests/testthat/test-stancode_snapshots.R")'

require(PDMPSamplersR, quietly = TRUE)

skip_if_no_stancode_snapshots <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not(
    identical(Sys.getenv("PDMPSAMPLERSR_EXHAUSTIVE_TESTS"), "true"),
    "Stan code snapshot tests not enabled (set PDMPSAMPLERSR_EXHAUSTIVE_TESTS=true)"
  )
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("rstan")
}

# Helper: generate the custom-family stancode (pre-rewrite) and extract
# the model block.  This is what the rewrite functions actually operate on.
pre_rewrite_model_block <- function(formula, data, family, ...) {
  sub_family <- PDMPSamplersR:::make_subsampled_family(family)
  sub_stanvars <- PDMPSamplersR:::make_pdmp_stanvars(
    formula, data, family, ...
  )
  code <- brms::stancode(formula, data = data, family = sub_family,
                          stanvars = sub_stanvars, ...)
  PDMPSamplersR:::extract_named_block(code, "model")
}

# Helper: generate ext_cpp stancode and extract the rewritten model block
rewritten_model_block <- function(formula, data, family, ...) {
  result <- brm_stancode(formula, data = data, family = family,
                          subsample_size = 10L, ...)
  PDMPSamplersR:::extract_named_block(result$ext_cpp, "model")
}

# -- Shared data ---------------------------------------------------------------

make_fe_data <- function() {
  set.seed(42)
  data.frame(y = rnorm(50), x = rnorm(50))
}

make_spline_data <- function() {
  set.seed(42)
  data.frame(y = rnorm(60), x = rnorm(60),
             z = seq(0, 2 * pi, length.out = 60))
}

make_gp_data <- function() {
  set.seed(42)
  data.frame(y = rnorm(50), x = seq(0, 4, length.out = 50))
}

make_re_data <- function() {
  set.seed(42)
  data.frame(y = rnorm(50), x = rnorm(50), group = rep(1:5, each = 10))
}

make_multi_group_data <- function() {
  set.seed(42)
  data.frame(y = rnorm(60), x = rnorm(60),
             g1 = rep(1:3, each = 20), g2 = rep(1:4, times = 15))
}

make_spline_re_data <- function() {
  set.seed(42)
  data.frame(y = rnorm(60), x = rnorm(60),
             z = seq(0, 2 * pi, length.out = 60),
             group = rep(1:6, each = 10))
}

# ==============================================================================
# Section 1: Low-level rewrite function snapshots
# ==============================================================================

test_that("rewrite snapshot: rewrite_mu_init on FE model", {
  skip_if_no_stancode_snapshots()
  block <- pre_rewrite_model_block(y ~ x, make_fe_data(), gaussian())
  result <- PDMPSamplersR:::rewrite_mu_init(block)
  expect_snapshot(cat(result))
})

test_that("rewrite snapshot: rewrite_fe_accumulation on FE model", {
  skip_if_no_stancode_snapshots()
  block <- pre_rewrite_model_block(y ~ x, make_fe_data(), gaussian())
  result <- PDMPSamplersR:::rewrite_fe_accumulation(block)
  expect_snapshot(cat(result))
})

test_that("rewrite snapshot: rewrite_spline_matrices on spline model", {
  skip_if_no_stancode_snapshots()
  block <- pre_rewrite_model_block(y ~ x + s(z, k = 5), make_spline_data(), gaussian())
  result <- PDMPSamplersR:::rewrite_spline_matrices(block)
  expect_snapshot(cat(result))
})

test_that("rewrite snapshot: rewrite_gp_indexing on GP model", {
  skip_if_no_stancode_snapshots()
  block <- pre_rewrite_model_block(y ~ gp(x), make_gp_data(), gaussian())
  result <- PDMPSamplersR:::rewrite_gp_indexing(block)
  expect_snapshot(cat(result))
})

test_that("rewrite snapshot: rewrite_re_loops on RE intercept model", {
  skip_if_no_stancode_snapshots()
  block <- pre_rewrite_model_block(y ~ x + (1 | group), make_re_data(), gaussian())
  result <- PDMPSamplersR:::rewrite_re_loops(block)
  expect_snapshot(cat(result))
})

test_that("rewrite snapshot: rewrite_re_loops on RE slope model", {
  skip_if_no_stancode_snapshots()
  block <- pre_rewrite_model_block(
    y ~ x + (1 + x | group), make_re_data(), gaussian()
  )
  result <- PDMPSamplersR:::rewrite_re_loops(block)
  expect_snapshot(cat(result))
})

test_that("rewrite snapshot: rewrite_re_loops on multi-group RE model", {
  skip_if_no_stancode_snapshots()
  block <- pre_rewrite_model_block(
    y ~ x + (1 | g1) + (1 | g2), make_multi_group_data(), gaussian()
  )
  result <- PDMPSamplersR:::rewrite_re_loops(block)
  expect_snapshot(cat(result))
})

# ==============================================================================
# Section 2: Full pipeline model block snapshots
# ==============================================================================

test_that("model block snapshot: gaussian y ~ x", {
  skip_if_no_stancode_snapshots()
  block <- rewritten_model_block(y ~ x, make_fe_data(), gaussian())
  expect_snapshot(cat(block))
})

test_that("model block snapshot: gaussian y ~ x + s(z)", {
  skip_if_no_stancode_snapshots()
  block <- rewritten_model_block(
    y ~ x + s(z, k = 5), make_spline_data(), gaussian()
  )
  expect_snapshot(cat(block))
})

test_that("model block snapshot: gaussian y ~ gp(x)", {
  skip_if_no_stancode_snapshots()
  block <- rewritten_model_block(y ~ gp(x), make_gp_data(), gaussian())
  expect_snapshot(cat(block))
})

test_that("model block snapshot: gaussian y ~ x + (1 | group)", {
  skip_if_no_stancode_snapshots()
  block <- rewritten_model_block(
    y ~ x + (1 | group), make_re_data(), gaussian()
  )
  expect_snapshot(cat(block))
})

test_that("model block snapshot: bernoulli y ~ x + (1 | group)", {
  skip_if_no_stancode_snapshots()
  set.seed(42)
  df <- data.frame(y = rbinom(50, 1, 0.5), x = rnorm(50),
                   group = rep(1:5, each = 10))
  block <- rewritten_model_block(y ~ x + (1 | group), df, brms::bernoulli())
  expect_snapshot(cat(block))
})

test_that("model block snapshot: gaussian y ~ x + (1 + x | group)", {
  skip_if_no_stancode_snapshots()
  block <- rewritten_model_block(
    y ~ x + (1 + x | group), make_re_data(), gaussian()
  )
  expect_snapshot(cat(block))
})

test_that("model block snapshot: gaussian y ~ x + (1 | g1) + (1 | g2)", {
  skip_if_no_stancode_snapshots()
  block <- rewritten_model_block(
    y ~ x + (1 | g1) + (1 | g2), make_multi_group_data(), gaussian()
  )
  expect_snapshot(cat(block))
})

test_that("model block snapshot: gaussian y ~ x + s(z) + (1 | group)", {
  skip_if_no_stancode_snapshots()
  block <- rewritten_model_block(
    y ~ x + s(z, k = 5) + (1 | group), make_spline_re_data(), gaussian()
  )
  expect_snapshot(cat(block))
})
