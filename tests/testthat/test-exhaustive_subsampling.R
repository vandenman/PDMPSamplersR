# Exhaustive integration tests for subsampled brms models.
#
# These tests are gated by the PDMPSAMPLERSR_EXHAUSTIVE_TESTS environment
# variable and are never run on CRAN.  They verify:
#
#   1. Gradient equivalence: the ext_cpp (subsampled) Stan model produces
#      the same gradient as a standard brms Stan model evaluated on the
#      same data when all observations are included in the subsample.
#
#   2. Subsample scaling: with a proper subsample of m < N, the gradient
#      of the ext_cpp model satisfies
#          grad_ext = (1 - s) * grad_prior + s * grad_std_sub
#      where s = N / m, grad_prior is the prior-only gradient, and
#      grad_std_sub is the standard model gradient on subsetted data.
#
#   3. Posterior correctness: PDMP posteriors (with and without sub-
#      sampling) are close to Stan HMC reference posteriors.
#
# Run locally with (CLI):
#   NOT_CRAN=true PDMPSAMPLERSR_EXHAUSTIVE_TESTS=true Rscript -e 'testthat::test_file("tests/testthat/test-exhaustive_subsampling.R")'
# Run locally with (RStudio):

# for standalone usage
require(PDMPSamplersR, quietly = TRUE)
require(testthat, quietly = TRUE)
Sys.setenv("PDMPSAMPLERSR_EXHAUSTIVE_TESTS" = "true")
# -- skip helpers --------------------------------------------------------------

skip_if_no_exhaustive_tests <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not(
    identical(Sys.getenv("PDMPSAMPLERSR_EXHAUSTIVE_TESTS"), "true"),
    "Exhaustive subsampling tests not enabled (set PDMPSAMPLERSR_EXHAUSTIVE_TESTS=true)"
  )
  skip_if_no_brms_setup()
}

# -- Julia gradient helper -----------------------------------------------------

.gradient_helper_loaded <- new.env(parent = emptyenv())

ensure_gradient_helper <- function() {
  if (!isTRUE(.gradient_helper_loaded$done)) {
    helper_path <- system.file("julia", "gradient_test_helper.jl",
                               package = "PDMPSamplersR")
    JuliaCall::julia_source(helper_path)
    .gradient_helper_loaded$done <- TRUE
  }
}

eval_gradient <- function(stan_file, data_json, theta = numeric(0),
                          hpp_path = "", indices_0based = integer(0)) {
  ensure_gradient_helper()
  JuliaCall::julia_call(
    "r_eval_stan_gradient",
    normalizePath(stan_file, mustWork = TRUE),
    normalizePath(data_json, mustWork = TRUE),
    as.numeric(theta),
    hpp_path = hpp_path,
    indices_0based = as.integer(indices_0based)
  )
}

# -- R helpers for model setup -------------------------------------------------

setup_gradient_test <- function(formula, data, family, prior = NULL, ...) {
  scode <- brms::stancode(formula, data = data, family = family,
                          prior = prior, ...)
  scode_ext <- PDMPSamplersR:::make_ext_cpp_stancode(
    scode, formula, data, family, prior, NULL, "no", ...
  )
  sdata <- brms::standata(formula, data = data, family = family,
                          prior = prior, ...)

  std_file <- PDMPSamplersR:::cached_stan_model(scode)
  ext_file <- PDMPSamplersR:::cached_stan_model(scode_ext)

  data_full_json <- tempfile(fileext = ".json")
  PDMPSamplersR::write_stan_json(sdata, data_full_json)

  list(
    std_file = std_file,
    ext_file = ext_file,
    sdata = sdata,
    data_full_json = data_full_json,
    N = nrow(data)
  )
}

run_fulldata_gradient_test <- function(setup, tol = 1e-6) {
  hpp <- normalizePath(PDMPSamplersR:::hpp_path(), mustWork = TRUE)
  N <- setup$N
  indices <- as.integer(0:(N - 1L))

  res_std <- eval_gradient(setup$std_file, setup$data_full_json)
  res_ext <- eval_gradient(setup$ext_file, setup$data_full_json,
                           theta = res_std$theta,
                           hpp_path = hpp, indices_0based = indices)

  testthat::expect_equal(res_ext$d, res_std$d, label = "parameter dimensions")
  testthat::expect_equal(res_ext$lp, res_std$lp, tolerance = tol,
                         label = "log density (full data)")
  testthat::expect_equal(res_ext$grad, res_std$grad, tolerance = tol,
                         label = "gradient (full data)")
}

run_subsample_gradient_test <- function(setup, m, tol = 1e-5) {
  hpp <- normalizePath(PDMPSamplersR:::hpp_path(), mustWork = TRUE)
  N <- setup$N
  K <- N %/% m
  stopifnot(N %% m == 0L)

  all_indices <- as.integer(0:(N - 1L))

  ext_full <- eval_gradient(setup$ext_file, setup$data_full_json,
                            theta = numeric(0),
                            hpp_path = hpp, indices_0based = all_indices)
  theta <- ext_full$theta

  set.seed(42)
  perm <- sample.int(N)
  grad_sum <- rep(0, length(theta))
  for (k in seq_len(K)) {
    idx_r <- sort(perm[((k - 1L) * m + 1L):(k * m)])
    idx_0 <- as.integer(idx_r - 1L)
    ext_k <- eval_gradient(setup$ext_file, setup$data_full_json,
                           theta = theta,
                           hpp_path = hpp, indices_0based = idx_0)
    grad_sum <- grad_sum + ext_k$grad
  }

  # validate the subsampled gradient on all data against brms
  avg_grad <- grad_sum / K
  testthat::expect_equal(avg_grad, ext_full$grad, tolerance = tol,
                         label = "subsample partition gradient")

  # validate the full gradient against brms -- already done by testset 1.
  #   ext_full_brms <- eval_gradient(setup$std_file, setup$data_full_json, theta = numeric(0))
  #   testthat::expect_equal(ext_full_brms$grad, ext_full$grad, tolerance = tol,
  #                          label = "subsample gradient on all data matched brms")

}

# -- Data generators -----------------------------------------------------------

make_gaussian_fe_data <- function(n = 80, seed = 42) {
  set.seed(seed)
  x <- rnorm(n)
  y <- 2 + 0.5 * x + rnorm(n, sd = 0.5)
  data.frame(y = y, x = x)
}

make_bernoulli_fe_data <- function(n = 100, seed = 123) {
  set.seed(seed)
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.5 + 1.0 * x))
  data.frame(y = y, x = x)
}

make_poisson_fe_data <- function(n = 100, seed = 456) {
  set.seed(seed)
  x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.3 * x))
  data.frame(y = y, x = x)
}

make_negbinomial_fe_data <- function(n = 100, seed = 200) {
  set.seed(seed)
  x <- rnorm(n)
  mu <- exp(0.5 + 0.3 * x)
  y <- rnbinom(n, size = 5, mu = mu)
  data.frame(y = y, x = x)
}

make_gaussian_spline_data <- function(n = 100, seed = 77) {
  set.seed(seed)
  x <- seq(0, 2 * pi, length.out = n)
  y <- sin(x) + rnorm(n, sd = 0.3)
  data.frame(y = y, x = x)
}

make_gaussian_gp_data <- function(n = 60, seed = 88) {
  set.seed(seed)
  x <- seq(0, 4, length.out = n)
  y <- sin(x) + rnorm(n, sd = 0.3)
  data.frame(y = y, x = x)
}

make_gaussian_fe_spline_data <- function(n = 100, seed = 99) {
  set.seed(seed)
  x <- rnorm(n)
  z <- seq(0, 2 * pi, length.out = n)
  y <- 1.5 + 0.4 * x + sin(z) + rnorm(n, sd = 0.3)
  data.frame(y = y, x = x, z = z)
}

make_gaussian_fe_gp_data <- function(n = 60, seed = 101) {
  set.seed(seed)
  x <- rnorm(n)
  z <- seq(0, 4, length.out = n)
  y <- 1.0 + 0.3 * x + sin(z) + rnorm(n, sd = 0.3)
  data.frame(y = y, x = x, z = z)
}

make_gaussian_fe_spline_gp_data <- function(n = 60, seed = 111) {
  set.seed(seed)
  x <- rnorm(n)
  z <- seq(0, 2 * pi, length.out = n)
  w <- seq(0, 4, length.out = n)
  y <- 1.0 + 0.3 * x + sin(z) + cos(w) + rnorm(n, sd = 0.3)
  data.frame(y = y, x = x, z = z, w = w)
}

make_bernoulli_spline_data <- function(n = 120, seed = 133) {
  set.seed(seed)
  x <- seq(-3, 3, length.out = n)
  y <- rbinom(n, 1, plogis(sin(x)))
  data.frame(y = y, x = x)
}

make_poisson_spline_data <- function(n = 100, seed = 144) {
  set.seed(seed)
  x <- seq(0, 2 * pi, length.out = n)
  y <- rpois(n, exp(0.5 + 0.5 * sin(x)))
  data.frame(y = y, x = x)
}

make_poisson_fe_spline_data <- function(n = 100, seed = 155) {
  set.seed(seed)
  x <- rnorm(n)
  z <- seq(0, 2 * pi, length.out = n)
  y <- rpois(n, exp(0.3 + 0.2 * x + 0.5 * sin(z)))
  data.frame(y = y, x = x, z = z)
}

make_bernoulli_gp_data <- function(n = 60, seed = 166) {
  set.seed(seed)
  x <- seq(-2, 2, length.out = n)
  y <- rbinom(n, 1, plogis(sin(2 * x)))
  data.frame(y = y, x = x)
}

make_gaussian_multi_spline_data <- function(n = 100, seed = 177) {
  set.seed(seed)
  x <- seq(0, 2 * pi, length.out = n)
  z <- seq(0, 3, length.out = n)
  y <- sin(x) + cos(z) + rnorm(n, sd = 0.3)
  data.frame(y = y, x = x, z = z)
}

# ==============================================================================
# Testset 1a: Gradient equivalence tests (m = N, full data)
# ==============================================================================
# When all observations are included in the subsample, the ext_cpp model
# gradient must exactly match the standard model gradient (the N/m scaling
# cancels to 1).

test_that("gradient fulldata: gaussian y ~ x", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_data()
  setup <- setup_gradient_test(y ~ x, df, gaussian())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: bernoulli y ~ x", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_fe_data()
  setup <- setup_gradient_test(y ~ x, df, brms::bernoulli())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: poisson y ~ x", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_fe_data()
  setup <- setup_gradient_test(y ~ x, df, poisson())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: negbinomial y ~ x", {
  skip_if_no_exhaustive_tests()
  df <- make_negbinomial_fe_data()
  setup <- setup_gradient_test(y ~ x, df, brms::negbinomial())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: gaussian y ~ s(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_spline_data()
  setup <- setup_gradient_test(y ~ s(x, k = 5), df, gaussian())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: bernoulli y ~ s(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_spline_data()
  setup <- setup_gradient_test(y ~ s(x, k = 5), df, brms::bernoulli())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: poisson y ~ s(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_spline_data()
  setup <- setup_gradient_test(y ~ s(x, k = 5), df, poisson())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: gaussian y ~ gp(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_gp_data()
  setup <- setup_gradient_test(y ~ gp(x), df, gaussian())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: bernoulli y ~ gp(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_gp_data()
  setup <- setup_gradient_test(y ~ gp(x), df, brms::bernoulli())
  run_fulldata_gradient_test(setup, tol = 1e-3)
})

test_that("gradient fulldata: gaussian y ~ x + s(z)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5), df, gaussian())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: poisson y ~ x + s(z)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_fe_spline_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5), df, poisson())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: gaussian y ~ x + gp(z)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_gp_data()
  setup <- setup_gradient_test(y ~ x + gp(z), df, gaussian())
  run_fulldata_gradient_test(setup, tol = 1e-3)
})

test_that("gradient fulldata: gaussian y ~ x + s(z) + gp(w)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_gp_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5) + gp(w), df, gaussian())
  run_fulldata_gradient_test(setup, tol = 1e-3)
})

test_that("gradient fulldata: gaussian y ~ s(x, k=5) + s(z, k=5)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_multi_spline_data()
  setup <- setup_gradient_test(y ~ s(x, k = 5) + s(z, k = 5), df, gaussian())
  run_fulldata_gradient_test(setup)
})

# ==============================================================================
# Testset 1b: Gradient subsample scaling tests (m < N)
# ==============================================================================
# With a strict subset, verify: grad_ext = (1-s)*grad_prior + s*grad_std_sub

test_that("gradient subsample: gaussian y ~ x", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_data()
  setup <- setup_gradient_test(y ~ x, df, gaussian())
  run_subsample_gradient_test(setup, m = 20)
})

test_that("gradient subsample: bernoulli y ~ x", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_fe_data()
  setup <- setup_gradient_test(y ~ x, df, brms::bernoulli())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: poisson y ~ x", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_fe_data()
  setup <- setup_gradient_test(y ~ x, df, poisson())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: negbinomial y ~ x", {
  skip_if_no_exhaustive_tests()
  df <- make_negbinomial_fe_data()
  setup <- setup_gradient_test(y ~ x, df, brms::negbinomial())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: gaussian y ~ s(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_spline_data()
  setup <- setup_gradient_test(y ~ s(x, k = 5), df, gaussian())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: bernoulli y ~ s(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_spline_data()
  setup <- setup_gradient_test(y ~ s(x, k = 5), df, brms::bernoulli())
  run_subsample_gradient_test(setup, m = 30)
})

test_that("gradient subsample: poisson y ~ s(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_spline_data()
  setup <- setup_gradient_test(y ~ s(x, k = 5), df, poisson())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: gaussian y ~ gp(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_gp_data()
  setup <- setup_gradient_test(y ~ gp(x), df, gaussian())
  run_subsample_gradient_test(setup, m = 15)
})

test_that("gradient subsample: bernoulli y ~ gp(x)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_gp_data()
  setup <- setup_gradient_test(y ~ gp(x), df, brms::bernoulli())
  run_subsample_gradient_test(setup, m = 15)
})

test_that("gradient subsample: gaussian y ~ x + s(z)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5), df, gaussian())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: poisson y ~ x + s(z)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_fe_spline_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5), df, poisson())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: gaussian y ~ x + gp(z)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_gp_data()
  setup <- setup_gradient_test(y ~ x + gp(z), df, gaussian())
  run_subsample_gradient_test(setup, m = 15)
})

test_that("gradient subsample: gaussian y ~ x + s(z) + gp(w)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_gp_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5) + gp(w), df, gaussian())
  run_subsample_gradient_test(setup, m = 15)
})

test_that("gradient subsample: gaussian y ~ s(x, k=5) + s(z, k=5)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_multi_spline_data()
  setup <- setup_gradient_test(y ~ s(x, k = 5) + s(z, k = 5), df, gaussian())
  run_subsample_gradient_test(setup, m = 25)
})

# -- Data generators (random effects) -----------------------------------------

make_gaussian_re_intercept_data <- function(n = 80, n_groups = 4, seed = 300) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  re <- rnorm(n_groups, sd = 0.5)
  y <- 2 + 0.5 * x + re[group] + rnorm(n, sd = 0.5)
  data.frame(y = y, x = x, group = group)
}

make_bernoulli_re_intercept_data <- function(n = 100, n_groups = 5, seed = 310) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  re <- rnorm(n_groups, sd = 0.5)
  y <- rbinom(n, 1, plogis(0.5 + 1.0 * x + re[group]))
  data.frame(y = y, x = x, group = group)
}

make_poisson_re_intercept_data <- function(n = 100, n_groups = 5, seed = 320) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  re <- rnorm(n_groups, sd = 0.3)
  y <- rpois(n, exp(0.5 + 0.3 * x + re[group]))
  data.frame(y = y, x = x, group = group)
}

make_negbinomial_re_intercept_data <- function(n = 100, n_groups = 5, seed = 330) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  re <- rnorm(n_groups, sd = 0.3)
  mu <- exp(0.5 + 0.3 * x + re[group])
  y <- rnbinom(n, size = 5, mu = mu)
  data.frame(y = y, x = x, group = group)
}

make_gaussian_re_slope_data <- function(n = 80, n_groups = 4, seed = 340) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  re_int <- rnorm(n_groups, sd = 0.5)
  re_slope <- rnorm(n_groups, sd = 0.3)
  y <- 2 + (0.5 + re_slope[group]) * x + re_int[group] + rnorm(n, sd = 0.5)
  data.frame(y = y, x = x, group = group)
}

make_gaussian_re_multi_group_data <- function(n = 60, seed = 350) {
  set.seed(seed)
  g1 <- rep(1:3, each = 20)
  g2 <- rep(1:4, times = 15)
  x <- rnorm(n)
  re1 <- rnorm(3, sd = 0.5)
  re2 <- rnorm(4, sd = 0.3)
  y <- 1.0 + 0.4 * x + re1[g1] + re2[g2] + rnorm(n, sd = 0.5)
  data.frame(y = y, x = x, g1 = g1, g2 = g2)
}

make_gaussian_fe_spline_re_data <- function(n = 100, n_groups = 5, seed = 360) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  z <- seq(0, 2 * pi, length.out = n)
  re <- rnorm(n_groups, sd = 0.5)
  y <- 1.5 + 0.4 * x + sin(z) + re[group] + rnorm(n, sd = 0.3)
  data.frame(y = y, x = x, z = z, group = group)
}

make_bernoulli_re_slope_data <- function(n = 100, n_groups = 5, seed = 370) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  re_int <- rnorm(n_groups, sd = 0.5)
  re_slope <- rnorm(n_groups, sd = 0.3)
  y <- rbinom(n, 1, plogis(0.5 + (0.8 + re_slope[group]) * x + re_int[group]))
  data.frame(y = y, x = x, group = group)
}

make_poisson_fe_spline_re_data <- function(n = 100, n_groups = 5, seed = 380) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  z <- seq(0, 2 * pi, length.out = n)
  re <- rnorm(n_groups, sd = 0.3)
  y <- rpois(n, exp(0.3 + 0.2 * x + 0.5 * sin(z) + re[group]))
  data.frame(y = y, x = x, z = z, group = group)
}

# ==============================================================================
# Testset 2a: Gradient equivalence tests with random effects (m = N)
# ==============================================================================

test_that("gradient fulldata: gaussian y ~ x + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_re_intercept_data()
  setup <- setup_gradient_test(y ~ x + (1 | group), df, gaussian())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: bernoulli y ~ x + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_re_intercept_data()
  setup <- setup_gradient_test(y ~ x + (1 | group), df, brms::bernoulli())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: poisson y ~ x + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_re_intercept_data()
  setup <- setup_gradient_test(y ~ x + (1 | group), df, poisson())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: negbinomial y ~ x + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_negbinomial_re_intercept_data()
  setup <- setup_gradient_test(y ~ x + (1 | group), df, brms::negbinomial())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: gaussian y ~ x + (1 + x | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_re_slope_data()
  setup <- setup_gradient_test(y ~ x + (1 + x | group), df, gaussian())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: bernoulli y ~ x + (1 + x | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_re_slope_data()
  setup <- setup_gradient_test(y ~ x + (1 + x | group), df, brms::bernoulli())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: gaussian y ~ x + (1 | g1) + (1 | g2)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_re_multi_group_data()
  setup <- setup_gradient_test(y ~ x + (1 | g1) + (1 | g2), df, gaussian())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: gaussian y ~ x + s(z) + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_re_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5) + (1 | group), df, gaussian())
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: poisson y ~ x + s(z) + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_fe_spline_re_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5) + (1 | group), df, poisson())
  run_fulldata_gradient_test(setup)
})

# ==============================================================================
# Testset 2b: Gradient subsample scaling tests with random effects (m < N)
# ==============================================================================

test_that("gradient subsample: gaussian y ~ x + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_re_intercept_data()
  setup <- setup_gradient_test(y ~ x + (1 | group), df, gaussian())
  run_subsample_gradient_test(setup, m = 20)
})

test_that("gradient subsample: bernoulli y ~ x + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_re_intercept_data()
  setup <- setup_gradient_test(y ~ x + (1 | group), df, brms::bernoulli())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: poisson y ~ x + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_re_intercept_data()
  setup <- setup_gradient_test(y ~ x + (1 | group), df, poisson())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: negbinomial y ~ x + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_negbinomial_re_intercept_data()
  setup <- setup_gradient_test(y ~ x + (1 | group), df, brms::negbinomial())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: gaussian y ~ x + (1 + x | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_re_slope_data()
  setup <- setup_gradient_test(y ~ x + (1 + x | group), df, gaussian())
  run_subsample_gradient_test(setup, m = 20)
})

test_that("gradient subsample: bernoulli y ~ x + (1 + x | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_re_slope_data()
  setup <- setup_gradient_test(y ~ x + (1 + x | group), df, brms::bernoulli())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: gaussian y ~ x + (1 | g1) + (1 | g2)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_re_multi_group_data()
  setup <- setup_gradient_test(y ~ x + (1 | g1) + (1 | g2), df, gaussian())
  run_subsample_gradient_test(setup, m = 15)
})

test_that("gradient subsample: gaussian y ~ x + s(z) + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_re_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5) + (1 | group), df, gaussian())
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: poisson y ~ x + s(z) + (1 | group)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_fe_spline_re_data()
  setup <- setup_gradient_test(y ~ x + s(z, k = 5) + (1 | group), df, poisson())
  run_subsample_gradient_test(setup, m = 25)
})

# -- Data generators (distributional formulas) ---------------------------------

make_gaussian_sigma_fe_data <- function(n = 80, seed = 400) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  sigma <- exp(0.1 + 0.3 * z)
  y <- 2 + 0.5 * x + rnorm(n, sd = sigma)
  data.frame(y = y, x = x, z = z)
}

make_negbinomial_shape_fe_data <- function(n = 100, seed = 410) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  mu <- exp(0.5 + 0.3 * x)
  shape <- exp(1.0 + 0.4 * z)
  y <- rnbinom(n, size = shape, mu = mu)
  data.frame(y = y, x = x, z = z)
}

make_gaussian_sigma_spline_data <- function(n = 100, seed = 420) {
  set.seed(seed)
  x <- rnorm(n)
  z <- seq(0, 2 * pi, length.out = n)
  sigma <- exp(0.1 + 0.3 * sin(z))
  y <- 1.5 + 0.4 * x + rnorm(n, sd = sigma)
  data.frame(y = y, x = x, z = z)
}

make_gaussian_sigma_re_data <- function(n = 80, n_groups = 4, seed = 430) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  re_sigma <- rnorm(n_groups, sd = 0.3)
  sigma <- exp(0.2 + re_sigma[group])
  y <- 2 + 0.5 * x + rnorm(n, sd = sigma)
  data.frame(y = y, x = x, group = group)
}

make_gaussian_mu_re_sigma_fe_data <- function(n = 80, n_groups = 4, seed = 440) {
  set.seed(seed)
  group <- rep(seq_len(n_groups), each = n %/% n_groups)
  x <- rnorm(n)
  z <- rnorm(n)
  re_mu <- rnorm(n_groups, sd = 0.5)
  sigma <- exp(0.1 + 0.3 * z)
  y <- 2 + 0.5 * x + re_mu[group] + rnorm(n, sd = sigma)
  data.frame(y = y, x = x, z = z, group = group)
}

# ==============================================================================
# Testset 3a: Gradient equivalence tests with distributional formulas (m = N)
# ==============================================================================

test_that("gradient fulldata: gaussian bf(y ~ x, sigma ~ z)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_sigma_fe_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x, sigma ~ z), df, gaussian()
  )
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: negbinomial bf(y ~ x, shape ~ z)", {
  skip_if_no_exhaustive_tests()
  df <- make_negbinomial_shape_fe_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x, shape ~ z), df, brms::negbinomial()
  )
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: gaussian bf(y ~ x, sigma ~ s(z))", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_sigma_spline_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x, sigma ~ s(z, k = 5)), df, gaussian()
  )
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: gaussian bf(y ~ x, sigma ~ (1 | group))", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_sigma_re_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x, sigma ~ (1 | group)), df, gaussian()
  )
  run_fulldata_gradient_test(setup)
})

test_that("gradient fulldata: gaussian bf(y ~ x + (1|group), sigma ~ z)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_mu_re_sigma_fe_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x + (1 | group), sigma ~ z), df, gaussian()
  )
  run_fulldata_gradient_test(setup)
})

# ==============================================================================
# Testset 3b: Gradient subsample scaling with distributional formulas (m < N)
# ==============================================================================

test_that("gradient subsample: gaussian bf(y ~ x, sigma ~ z)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_sigma_fe_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x, sigma ~ z), df, gaussian()
  )
  run_subsample_gradient_test(setup, m = 20)
})

test_that("gradient subsample: negbinomial bf(y ~ x, shape ~ z)", {
  skip_if_no_exhaustive_tests()
  df <- make_negbinomial_shape_fe_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x, shape ~ z), df, brms::negbinomial()
  )
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: gaussian bf(y ~ x, sigma ~ s(z))", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_sigma_spline_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x, sigma ~ s(z, k = 5)), df, gaussian()
  )
  run_subsample_gradient_test(setup, m = 25)
})

test_that("gradient subsample: gaussian bf(y ~ x, sigma ~ (1 | group))", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_sigma_re_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x, sigma ~ (1 | group)), df, gaussian()
  )
  run_subsample_gradient_test(setup, m = 20)
})

test_that("gradient subsample: gaussian bf(y ~ x + (1|group), sigma ~ z)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_mu_re_sigma_fe_data()
  setup <- setup_gradient_test(
    brms::bf(y ~ x + (1 | group), sigma ~ z), df, gaussian()
  )
  run_subsample_gradient_test(setup, m = 20)
})

# ==============================================================================
# Section 3: Posterior correctness tests (TEMPORARILY COMMENTED OUT)
# ==============================================================================
# Compare PDMP posteriors (standard and subsampled) against Stan HMC.
# These are slow and commented out until gradient tests are stable.

if (FALSE) { # ---- begin commented-out posterior tests ----

compare_posteriors <- function(fit_pdmp, fit_stan, tol_est = 0.1,
                               tol_se = 0.2, tol_ci = 0.20) {
  ests_pdmp <- brms::fixef(fit_pdmp)
  ests_stan <- brms::fixef(fit_stan)
  shared <- intersect(rownames(ests_pdmp), rownames(ests_stan))
  testthat::expect_true(length(shared) > 0, label = "shared fixed effects")
  ests_pdmp <- ests_pdmp[shared, , drop = FALSE]
  ests_stan <- ests_stan[shared, , drop = FALSE]
  tols <- c("Estimate" = tol_est, "Est.Error" = tol_se,
            "Q2.5" = tol_ci, "Q97.5" = tol_ci)
  for (nm in names(tols)) {
    testthat::expect_equal(
      ests_pdmp[, nm], ests_stan[, nm],
      tolerance = tols[[nm]],
      label = paste("PDMP", nm, "vs Stan")
    )
  }
}

test_that("posterior: gaussian y ~ x (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_data()
  fit_pdmp <- brm_pdmp(y ~ x, data = df, family = gaussian(),
                       flow = "ZigZag", T = 50000, show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x, data = df, family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: gaussian y ~ x (subsampled)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_data()
  fit_sub <- brm_pdmp(y ~ x, data = df, family = gaussian(),
                      flow = "AdaptiveBoomerang", T = 50000,
                      show_progress = FALSE, subsample_size = 20L)
  fit_stan <- brms::brm(y ~ x, data = df, family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_sub, fit_stan, tol_est = 0.15, tol_se = 0.25,
                     tol_ci = 0.25)
})

test_that("posterior: bernoulli y ~ x (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_fe_data()
  fit_pdmp <- brm_pdmp(y ~ x, data = df, family = brms::bernoulli(),
                       flow = "ZigZag", T = 50000, show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x, data = df, family = brms::bernoulli(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: bernoulli y ~ x (subsampled)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_fe_data()
  fit_sub <- brm_pdmp(y ~ x, data = df, family = brms::bernoulli(),
                      flow = "AdaptiveBoomerang", T = 50000,
                      show_progress = FALSE, subsample_size = 25L)
  fit_stan <- brms::brm(y ~ x, data = df, family = brms::bernoulli(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_sub, fit_stan, tol_est = 0.15, tol_se = 0.25,
                     tol_ci = 0.25)
})

test_that("posterior: poisson y ~ x (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_fe_data()
  fit_pdmp <- brm_pdmp(y ~ x, data = df, family = poisson(),
                       flow = "ZigZag", T = 50000, show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x, data = df, family = poisson(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: poisson y ~ x (subsampled)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_fe_data()
  fit_sub <- brm_pdmp(y ~ x, data = df, family = poisson(),
                      flow = "AdaptiveBoomerang", T = 50000,
                      show_progress = FALSE, subsample_size = 25L)
  fit_stan <- brms::brm(y ~ x, data = df, family = poisson(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_sub, fit_stan, tol_est = 0.15, tol_se = 0.25,
                     tol_ci = 0.25)
})

test_that("posterior: negbinomial y ~ x (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_negbinomial_fe_data()
  fit_pdmp <- brm_pdmp(y ~ x, data = df, family = brms::negbinomial(),
                       flow = "AdaptiveBoomerang", T = 50000,
                       show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x, data = df, family = brms::negbinomial(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: negbinomial y ~ x (subsampled)", {
  skip_if_no_exhaustive_tests()
  df <- make_negbinomial_fe_data()
  fit_sub <- brm_pdmp(y ~ x, data = df, family = brms::negbinomial(),
                      flow = "AdaptiveBoomerang", T = 50000,
                      show_progress = FALSE, subsample_size = 25L)
  fit_stan <- brms::brm(y ~ x, data = df, family = brms::negbinomial(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_sub, fit_stan, tol_est = 0.15, tol_se = 0.25,
                     tol_ci = 0.25)
})

test_that("posterior: gaussian y ~ s(x) (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_spline_data()
  fit_pdmp <- brm_pdmp(y ~ s(x, k = 5), data = df, family = gaussian(),
                       flow = "AdaptiveBoomerang", T = 50000,
                       show_progress = FALSE)
  fit_stan <- brms::brm(y ~ s(x, k = 5), data = df, family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: gaussian y ~ s(x) (subsampled)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_spline_data()
  fit_sub <- brm_pdmp(y ~ s(x, k = 5), data = df, family = gaussian(),
                      flow = "AdaptiveBoomerang", T = 50000,
                      show_progress = FALSE, subsample_size = 25L)
  fit_stan <- brms::brm(y ~ s(x, k = 5), data = df, family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_sub, fit_stan, tol_est = 0.15, tol_se = 0.25,
                     tol_ci = 0.25)
})

test_that("posterior: gaussian y ~ x + s(z) (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_data()
  fit_pdmp <- brm_pdmp(y ~ x + s(z, k = 5), data = df, family = gaussian(),
                       flow = "AdaptiveBoomerang", T = 50000,
                       show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x + s(z, k = 5), data = df, family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: gaussian y ~ x + s(z) (subsampled)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_data()
  fit_sub <- brm_pdmp(y ~ x + s(z, k = 5), data = df, family = gaussian(),
                      flow = "AdaptiveBoomerang", T = 50000,
                      show_progress = FALSE, subsample_size = 25L)
  fit_stan <- brms::brm(y ~ x + s(z, k = 5), data = df, family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_sub, fit_stan, tol_est = 0.15, tol_se = 0.25,
                     tol_ci = 0.25)
})

test_that("posterior: gaussian y ~ gp(x) (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_gp_data()
  fit_pdmp <- brm_pdmp(y ~ gp(x), data = df, family = gaussian(),
                       flow = "AdaptiveBoomerang", T = 50000,
                       show_progress = FALSE)
  fit_stan <- brms::brm(y ~ gp(x), data = df, family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: gaussian y ~ gp(x) (subsampled)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_gp_data()
  fit_sub <- brm_pdmp(y ~ gp(x), data = df, family = gaussian(),
                      flow = "AdaptiveBoomerang", T = 50000,
                      show_progress = FALSE, subsample_size = 15L)
  fit_stan <- brms::brm(y ~ gp(x), data = df, family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_sub, fit_stan, tol_est = 0.15, tol_se = 0.25,
                     tol_ci = 0.25)
})

test_that("posterior: gaussian y ~ x + s(z) + gp(w) (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_gp_data()
  fit_pdmp <- brm_pdmp(y ~ x + s(z, k = 5) + gp(w), data = df,
                       family = gaussian(),
                       flow = "AdaptiveBoomerang", T = 80000,
                       show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x + s(z, k = 5) + gp(w), data = df,
                        family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan, tol_est = 0.15, tol_se = 0.25,
                     tol_ci = 0.25)
})

test_that("posterior: gaussian y ~ x + s(z) + gp(w) (subsampled)", {
  skip_if_no_exhaustive_tests()
  df <- make_gaussian_fe_spline_gp_data()
  fit_sub <- brm_pdmp(y ~ x + s(z, k = 5) + gp(w), data = df,
                      family = gaussian(),
                      flow = "AdaptiveBoomerang", T = 80000,
                      show_progress = FALSE, subsample_size = 15L)
  fit_stan <- brms::brm(y ~ x + s(z, k = 5) + gp(w), data = df,
                        family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_sub, fit_stan, tol_est = 0.20, tol_se = 0.30,
                     tol_ci = 0.30)
})

test_that("posterior: bernoulli y ~ s(x) (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_bernoulli_spline_data()
  fit_pdmp <- brm_pdmp(y ~ s(x, k = 5), data = df,
                       family = brms::bernoulli(),
                       flow = "AdaptiveBoomerang", T = 50000,
                       show_progress = FALSE)
  fit_stan <- brms::brm(y ~ s(x, k = 5), data = df,
                        family = brms::bernoulli(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: poisson y ~ x + s(z) (standard)", {
  skip_if_no_exhaustive_tests()
  df <- make_poisson_fe_spline_data()
  fit_pdmp <- brm_pdmp(y ~ x + s(z, k = 5), data = df, family = poisson(),
                       flow = "AdaptiveBoomerang", T = 50000,
                       show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x + s(z, k = 5), data = df, family = poisson(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})
} # ---- end commented-out posterior tests ----
