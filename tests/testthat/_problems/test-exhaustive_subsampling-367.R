# Extracted from test-exhaustive_subsampling.R:367

# prequel ----------------------------------------------------------------------
require(PDMPSamplersR)
skip_if_no_exhaustive_tests <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not(
    identical(Sys.getenv("PDMPSAMPLERSR_EXHAUSTIVE_TESTS"), "true"),
    "Exhaustive subsampling tests not enabled (set PDMPSAMPLERSR_EXHAUSTIVE_TESTS=true)"
  )
  skip_if_no_brms_setup()
}
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
  s <- N / m
  indices_r <- sort(sample.int(N, m))
  indices_0 <- as.integer(indices_r - 1L)

  sdata_sub <- PDMPSamplersR:::subset_standata(setup$sdata, indices_r)
  sdata_prior <- PDMPSamplersR:::make_prior_standata(setup$sdata)

  data_sub_json <- tempfile(fileext = ".json")
  data_prior_json <- tempfile(fileext = ".json")
  PDMPSamplersR::write_stan_json(sdata_sub, data_sub_json)
  PDMPSamplersR::write_stan_json(sdata_prior, data_prior_json)

  res_std_sub <- eval_gradient(setup$std_file, data_sub_json)
  theta <- res_std_sub$theta

  res_ext <- eval_gradient(setup$ext_file, setup$data_full_json,
                           theta = theta,
                           hpp_path = hpp, indices_0based = indices_0)
  res_prior <- eval_gradient(setup$std_file, data_prior_json, theta = theta)

  expected_grad <- (1 - s) * res_prior$grad + s * res_std_sub$grad
  testthat::expect_equal(res_ext$grad, expected_grad, tolerance = tol,
                         label = "subsample scaling gradient")
}
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

# test -------------------------------------------------------------------------
skip_if_no_exhaustive_tests()
df <- make_poisson_fe_data()
setup <- setup_gradient_test(y ~ x, df, poisson())
run_subsample_gradient_test(setup, m = 25)
