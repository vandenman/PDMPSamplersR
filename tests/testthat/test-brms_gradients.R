# Gradient equivalence tests for subsampled brms models.
#
# These tests verify:
#   1. Gradient equivalence: the ext_cpp (subsampled) Stan model produces
#      the same gradient as a standard brms Stan model on full data.
#   2. Subsample scaling: with m < N, the gradient partitions correctly
#      so that averaging over non-overlapping subsamples recovers the
#      full-data gradient.
#
# Requires Julia + BridgeStan (~22s total).  Always skipped on CRAN.

require(PDMPSamplersR, quietly = TRUE)
require(testthat, quietly = TRUE)

# -- skip helper ---------------------------------------------------------------

skip_if_no_gradient_tests <- function() {
  testthat::skip_on_cran()
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

make_student_fe_data <- function(n = 50, seed = 200) {
  set.seed(seed)
  x <- rnorm(n)
  y <- brms::rstudent_t(n, df = 5, mu = 1 + 0.5 * x, sigma = 1)
  data.frame(y = y, x = x)
}

make_gamma_fe_data <- function(n = 50, seed = 210) {
  set.seed(seed)
  x <- rnorm(n)
  mu <- exp(0.5 + 0.3 * x)
  shape <- 2
  y <- rgamma(n, shape = shape, rate = shape / mu)
  data.frame(y = y, x = x)
}

make_beta_fe_data <- function(n = 50, seed = 220) {
  set.seed(seed)
  x <- rnorm(n)
  mu <- plogis(0.2 + 0.5 * x)
  phi <- 10
  y <- rbeta(n, shape1 = mu * phi, shape2 = (1 - mu) * phi)
  data.frame(y = y, x = x)
}

make_student_sigma_fe_data <- function(n = 50, seed = 230) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  sigma <- exp(0.5 + 0.3 * z)
  y <- brms::rstudent_t(n, df = 5, mu = 1 + 0.5 * x, sigma = sigma)
  data.frame(y = y, x = x, z = z)
}

make_gamma_shape_fe_data <- function(n = 50, seed = 240) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  mu <- exp(0.5 + 0.3 * x)
  shape <- exp(0.7 + 0.2 * z)
  y <- rgamma(n, shape = shape, rate = shape / mu)
  data.frame(y = y, x = x, z = z)
}

# ==============================================================================
# Gradient test specs (table-driven)
# ==============================================================================
# Each spec defines a gradient test case. The fulldata and subsample loops
# iterate over these specs to generate test_that blocks automatically.
#
# Fields:
#   label    - test description (shown in testthat output)
#   data_fn  - function() returning a data.frame
#   formula  - brms formula
#   family   - brms family object
#   tol      - (optional) tolerance for fulldata test (default 1e-6)
#   m        - subsample size for subsample test (N must be divisible by m)

gradient_specs_fe <- list(
  list(label = "gaussian y ~ x",
       data_fn = make_gaussian_fe_data, formula = y ~ x,
       family = gaussian(), m = 20),
  list(label = "bernoulli y ~ x",
       data_fn = make_bernoulli_fe_data, formula = y ~ x,
       family = brms::bernoulli(), m = 25),
  list(label = "poisson y ~ x",
       data_fn = make_poisson_fe_data, formula = y ~ x,
       family = poisson(), m = 25),
  list(label = "negbinomial y ~ x",
       data_fn = make_negbinomial_fe_data, formula = y ~ x,
       family = brms::negbinomial(), m = 25),
  list(label = "gaussian y ~ s(x)",
       data_fn = make_gaussian_spline_data, formula = y ~ s(x, k = 5),
       family = gaussian(), m = 25),
  list(label = "bernoulli y ~ s(x)",
       data_fn = make_bernoulli_spline_data, formula = y ~ s(x, k = 5),
       family = brms::bernoulli(), m = 30),
  list(label = "poisson y ~ s(x)",
       data_fn = make_poisson_spline_data, formula = y ~ s(x, k = 5),
       family = poisson(), m = 25),
  list(label = "gaussian y ~ gp(x)",
       data_fn = make_gaussian_gp_data, formula = y ~ gp(x),
       family = gaussian(), tol = 1e-3, tol_sub = 1e-3, m = 15),
  list(label = "bernoulli y ~ gp(x)",
       data_fn = make_bernoulli_gp_data, formula = y ~ gp(x),
       family = brms::bernoulli(), tol = 1e-3, tol_sub = 1e-3, m = 15),
  list(label = "gaussian y ~ x + s(z)",
       data_fn = make_gaussian_fe_spline_data, formula = y ~ x + s(z, k = 5),
       family = gaussian(), m = 25),
  list(label = "poisson y ~ x + s(z)",
       data_fn = make_poisson_fe_spline_data, formula = y ~ x + s(z, k = 5),
       family = poisson(), m = 25),
  list(label = "gaussian y ~ x + gp(z)",
       data_fn = make_gaussian_fe_gp_data, formula = y ~ x + gp(z),
       family = gaussian(), tol = 1e-3, tol_sub = 1e-3, m = 15),
  list(label = "gaussian y ~ x + s(z) + gp(w)",
       data_fn = make_gaussian_fe_spline_gp_data,
       formula = y ~ x + s(z, k = 5) + gp(w),
       family = gaussian(), tol = 1e-3, tol_sub = 1e-3, m = 15),
  list(label = "gaussian y ~ s(x, k=5) + s(z, k=5)",
       data_fn = make_gaussian_multi_spline_data,
       formula = y ~ s(x, k = 5) + s(z, k = 5),
       family = gaussian(), m = 25),
  list(label = "student y ~ x",
       data_fn = make_student_fe_data, formula = y ~ x,
       family = brms::student(), m = 25),
  list(label = "Gamma y ~ x",
       data_fn = make_gamma_fe_data, formula = y ~ x,
       family = brms::brmsfamily("Gamma"), m = 25),
  list(label = "Beta y ~ x",
       data_fn = make_beta_fe_data, formula = y ~ x,
       family = brms::brmsfamily("Beta"), m = 25)
)

for (spec in gradient_specs_fe) {
  test_that(paste0("gradient fulldata: ", spec$label), {
    skip_if_no_gradient_tests()
    df <- spec$data_fn()
    setup <- setup_gradient_test(spec$formula, df, spec$family)
    run_fulldata_gradient_test(setup, tol = spec$tol %||% 1e-6)
  })
}

for (spec in gradient_specs_fe) {
  test_that(paste0("gradient subsample (m=", spec$m, "): ", spec$label), {
    skip_if_no_gradient_tests()
    df <- spec$data_fn()
    setup <- setup_gradient_test(spec$formula, df, spec$family)
    run_subsample_gradient_test(setup, m = spec$m, tol = spec$tol_sub %||% 1e-5)
  })
}

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

gradient_specs_re <- list(
  list(label = "gaussian y ~ x + (1 | group)",
       data_fn = make_gaussian_re_intercept_data,
       formula = y ~ x + (1 | group), family = gaussian(), m = 20),
  list(label = "bernoulli y ~ x + (1 | group)",
       data_fn = make_bernoulli_re_intercept_data,
       formula = y ~ x + (1 | group), family = brms::bernoulli(), m = 25),
  list(label = "poisson y ~ x + (1 | group)",
       data_fn = make_poisson_re_intercept_data,
       formula = y ~ x + (1 | group), family = poisson(), m = 25),
  list(label = "negbinomial y ~ x + (1 | group)",
       data_fn = make_negbinomial_re_intercept_data,
       formula = y ~ x + (1 | group), family = brms::negbinomial(), m = 25),
  list(label = "gaussian y ~ x + (1 + x | group)",
       data_fn = make_gaussian_re_slope_data,
       formula = y ~ x + (1 + x | group), family = gaussian(), m = 20),
  list(label = "bernoulli y ~ x + (1 + x | group)",
       data_fn = make_bernoulli_re_slope_data,
       formula = y ~ x + (1 + x | group), family = brms::bernoulli(), m = 25),
  list(label = "gaussian y ~ x + (1 | g1) + (1 | g2)",
       data_fn = make_gaussian_re_multi_group_data,
       formula = y ~ x + (1 | g1) + (1 | g2), family = gaussian(), m = 15),
  list(label = "gaussian y ~ x + s(z) + (1 | group)",
       data_fn = make_gaussian_fe_spline_re_data,
       formula = y ~ x + s(z, k = 5) + (1 | group), family = gaussian(), m = 25),
  list(label = "poisson y ~ x + s(z) + (1 | group)",
       data_fn = make_poisson_fe_spline_re_data,
       formula = y ~ x + s(z, k = 5) + (1 | group), family = poisson(), m = 25)
)

for (spec in gradient_specs_re) {
  test_that(paste0("gradient fulldata: ", spec$label), {
    skip_if_no_gradient_tests()
    df <- spec$data_fn()
    setup <- setup_gradient_test(spec$formula, df, spec$family)
    run_fulldata_gradient_test(setup, tol = spec$tol %||% 1e-6)
  })
}

for (spec in gradient_specs_re) {
  test_that(paste0("gradient subsample (m=", spec$m, "): ", spec$label), {
    skip_if_no_gradient_tests()
    df <- spec$data_fn()
    setup <- setup_gradient_test(spec$formula, df, spec$family)
    run_subsample_gradient_test(setup, m = spec$m, tol = spec$tol_sub %||% 1e-5)
  })
}

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

gradient_specs_dpar <- list(
  list(label = "gaussian bf(y ~ x, sigma ~ z)",
       data_fn = make_gaussian_sigma_fe_data,
       formula = brms::bf(y ~ x, sigma ~ z), family = gaussian(), m = 20),
  list(label = "negbinomial bf(y ~ x, shape ~ z)",
       data_fn = make_negbinomial_shape_fe_data,
       formula = brms::bf(y ~ x, shape ~ z), family = brms::negbinomial(), m = 25),
  list(label = "gaussian bf(y ~ x, sigma ~ s(z))",
       data_fn = make_gaussian_sigma_spline_data,
       formula = brms::bf(y ~ x, sigma ~ s(z, k = 5)), family = gaussian(), m = 25),
  list(label = "gaussian bf(y ~ x, sigma ~ (1 | group))",
       data_fn = make_gaussian_sigma_re_data,
       formula = brms::bf(y ~ x, sigma ~ (1 | group)), family = gaussian(), m = 20),
  list(label = "gaussian bf(y ~ x + (1|group), sigma ~ z)",
       data_fn = make_gaussian_mu_re_sigma_fe_data,
       formula = brms::bf(y ~ x + (1 | group), sigma ~ z), family = gaussian(), m = 20),
  list(label = "student bf(y ~ x, sigma ~ z)",
       data_fn = make_student_sigma_fe_data,
       formula = brms::bf(y ~ x, sigma ~ z), family = brms::student(), m = 25),
  list(label = "Gamma bf(y ~ x, shape ~ z)",
       data_fn = make_gamma_shape_fe_data,
       formula = brms::bf(y ~ x, shape ~ z), family = brms::brmsfamily("Gamma"), m = 25)
)

for (spec in gradient_specs_dpar) {
  test_that(paste0("gradient fulldata: ", spec$label), {
    skip_if_no_gradient_tests()
    df <- spec$data_fn()
    setup <- setup_gradient_test(spec$formula, df, spec$family)
    run_fulldata_gradient_test(setup, tol = spec$tol %||% 1e-6)
  })
}

for (spec in gradient_specs_dpar) {
  test_that(paste0("gradient subsample (m=", spec$m, "): ", spec$label), {
    skip_if_no_gradient_tests()
    df <- spec$data_fn()
    setup <- setup_gradient_test(spec$formula, df, spec$family)
    run_subsample_gradient_test(setup, m = spec$m, tol = spec$tol_sub %||% 1e-5)
  })
}

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
  skip_if_no_gradient_tests()
  df <- make_gaussian_fe_data()
  fit_pdmp <- brm_pdmp(y ~ x, data = df, family = gaussian(),
                       flow = "ZigZag", T = 50000, show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x, data = df, family = gaussian(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: gaussian y ~ x (subsampled)", {
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
  df <- make_bernoulli_fe_data()
  fit_pdmp <- brm_pdmp(y ~ x, data = df, family = brms::bernoulli(),
                       flow = "ZigZag", T = 50000, show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x, data = df, family = brms::bernoulli(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: bernoulli y ~ x (subsampled)", {
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
  df <- make_poisson_fe_data()
  fit_pdmp <- brm_pdmp(y ~ x, data = df, family = poisson(),
                       flow = "ZigZag", T = 50000, show_progress = FALSE)
  fit_stan <- brms::brm(y ~ x, data = df, family = poisson(),
                        chains = 2, iter = 2000, warmup = 1000,
                        silent = 2, refresh = 0)
  compare_posteriors(fit_pdmp, fit_stan)
})

test_that("posterior: poisson y ~ x (subsampled)", {
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
  skip_if_no_gradient_tests()
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
