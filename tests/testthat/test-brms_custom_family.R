# Tests for custom-family subsampling: family construction, Stan function
# generation, brace-aware helpers, and Phase 1 rewrites.

skip_if_no_brms <- function() {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("rstan")
}

# -- Family helpers ------------------------------------------------------------

test_that("subsampled_family_name maps baked-in links", {
  expect_equal(
    PDMPSamplersR:::subsampled_family_name("bernoulli", "logit"),
    "subsampled_bernoulli_logit"
  )
  expect_equal(
    PDMPSamplersR:::subsampled_family_name("poisson", "log"),
    "subsampled_poisson_log"
  )
  expect_equal(
    PDMPSamplersR:::subsampled_family_name("negbinomial", "log"),
    "subsampled_negbinomial_log"
  )
})

test_that("subsampled_family_name handles non-baked-in links", {
  expect_equal(
    PDMPSamplersR:::subsampled_family_name("gaussian", "identity"),
    "subsampled_gaussian"
  )
})

test_that("get_subsampled_family_links sets mu to identity", {
  skip_if_no_brms()

  fam <- brms::bernoulli()
  links <- PDMPSamplersR:::get_subsampled_family_links(fam)
  expect_equal(links, "identity")
})

test_that("get_subsampled_family_links preserves non-mu links", {
  skip_if_no_brms()

  fam <- gaussian()
  links <- PDMPSamplersR:::get_subsampled_family_links(fam)
  expect_equal(links[1], "identity")
  expect_equal(links[2], "log")
})

test_that("make_subsampled_family returns custom_family object", {
  skip_if_no_brms()

  sub_fam <- PDMPSamplersR:::make_subsampled_family(brms::bernoulli())
  expect_equal(sub_fam$name, "subsampled_bernoulli_logit")
  expect_equal(sub_fam$dpars, "mu")
  expect_false(sub_fam$loop)
  expect_equal(sub_fam$type, "int")
})

test_that("make_subsampled_family accepts student family", {
  skip_if_no_brms()

  sub_fam <- PDMPSamplersR:::make_subsampled_family(brms::student())
  expect_equal(sub_fam$name, "subsampled_student")
  expect_equal(sub_fam$dpars, c("mu", "sigma", "nu"))
  expect_false(sub_fam$loop)
  expect_equal(sub_fam$type, "real")
})

test_that("make_subsampled_family accepts Gamma family", {
  skip_if_no_brms()

  sub_fam <- PDMPSamplersR:::make_subsampled_family(brms::brmsfamily("Gamma"))
  expect_equal(sub_fam$name, "subsampled_gamma")
  expect_equal(sub_fam$dpars, c("mu", "shape"))
  expect_false(sub_fam$loop)
  expect_equal(sub_fam$type, "real")
})

test_that("make_subsampled_family accepts Beta family", {
  skip_if_no_brms()

  sub_fam <- PDMPSamplersR:::make_subsampled_family(brms::brmsfamily("Beta"))
  expect_equal(sub_fam$name, "subsampled_beta")
  expect_equal(sub_fam$dpars, c("mu", "phi"))
  expect_false(sub_fam$loop)
  expect_equal(sub_fam$type, "real")
})

test_that("get_subsampled_family_links preserves log link for Gamma mu", {
  skip_if_no_brms()

  fam <- brms::brmsfamily("Gamma")
  links <- PDMPSamplersR:::get_subsampled_family_links(fam)
  expect_equal(links[1], "log")
  expect_equal(links[2], "log")
})

test_that("get_subsampled_family_links preserves logit link for Beta mu", {
  skip_if_no_brms()

  fam <- brms::brmsfamily("Beta")
  links <- PDMPSamplersR:::get_subsampled_family_links(fam)
  expect_equal(links[1], "logit")
  expect_equal(links[2], "log")
})

test_that("get_subsampled_family_links uses identity for student mu", {
  skip_if_no_brms()

  fam <- brms::student()
  links <- PDMPSamplersR:::get_subsampled_family_links(fam)
  expect_equal(links[1], "identity")
  expect_equal(links[2], "log")
  expect_equal(links[3], "logm1")
})

test_that("get_subsampled_family_bounds includes nu and phi", {
  expect_equal(
    PDMPSamplersR:::get_subsampled_family_bounds(c("mu", "sigma", "nu")),
    c(NA, 0, 1)
  )
  expect_equal(
    PDMPSamplersR:::get_subsampled_family_bounds(c("mu", "phi")),
    c(NA, 0)
  )
})

# -- Stan function generation --------------------------------------------------

test_that("subsampled bernoulli_logit function uses local mu indexing", {
  codegen_info <- list(dpar_shapes = list())
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    brms::bernoulli(), codegen_info
  )
  expect_match(code, "subsampled_bernoulli_logit_lpmf")
  expect_match(code, "pdmp_get_subsample_size")
  expect_match(code, "pdmp_get_subsample_index")
  expect_match(code, "bernoulli_logit_lpmf")
  expect_match(code, "mu\\[i\\]")
  expect_no_match(code, "mu\\[idx\\]")
})

test_that("subsampled poisson_log function has correct structure", {
  codegen_info <- list(dpar_shapes = list())
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    poisson(), codegen_info
  )
  expect_match(code, "subsampled_poisson_log_lpmf")
  expect_match(code, "poisson_log_lpmf")
  expect_match(code, "mu\\[i\\]")
  expect_match(code, "y\\[idx\\]")
})

test_that("subsampled gaussian function includes scalar sigma", {
  codegen_info <- list(dpar_shapes = list(sigma = "scalar"))
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    gaussian(), codegen_info
  )
  expect_match(code, "subsampled_gaussian_lpdf")
  expect_match(code, "normal_lpdf")
  expect_match(code, "real sigma")
  expect_no_match(code, "vector sigma")
  expect_match(code, "mu\\[i\\], sigma")
})

test_that("subsampled gaussian function includes vector sigma", {
  codegen_info <- list(dpar_shapes = list(sigma = "vector"))
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    gaussian(), codegen_info
  )
  expect_match(code, "vector sigma")
  expect_match(code, "sigma\\[i\\]")
  expect_no_match(code, "sigma\\[idx\\]")
})

test_that("subsampled negbinomial_log function includes scalar shape", {
  codegen_info <- list(dpar_shapes = list(shape = "scalar"))
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    brms::negbinomial(), codegen_info
  )
  expect_match(code, "subsampled_negbinomial_log_lpmf")
  expect_match(code, "neg_binomial_2_log_lpmf")
  expect_match(code, "real shape")
  expect_match(code, "mu\\[i\\], shape\\b")
})

test_that("subsampled student function has correct call pattern", {
  codegen_info <- list(dpar_shapes = list(sigma = "scalar", nu = "scalar"))
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    brms::student(), codegen_info
  )
  expect_match(code, "subsampled_student_lpdf")
  expect_match(code, "student_t_lpdf")
  expect_match(code, "real sigma")
  expect_match(code, "real nu")
  expect_match(code, "nu, mu\\[i\\], sigma")
})

test_that("subsampled student function handles vector sigma", {
  codegen_info <- list(dpar_shapes = list(sigma = "vector", nu = "scalar"))
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    brms::student(), codegen_info
  )
  expect_match(code, "vector sigma")
  expect_match(code, "nu, mu\\[i\\], sigma\\[i\\]")
})

test_that("subsampled gamma function has correct call pattern", {
  codegen_info <- list(dpar_shapes = list(shape = "scalar"))
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    brms::brmsfamily("Gamma"), codegen_info
  )
  expect_match(code, "subsampled_gamma_lpdf")
  expect_match(code, "gamma_lpdf")
  expect_match(code, "real shape")
  expect_match(code, "shape, shape / mu\\[i\\]")
})

test_that("subsampled gamma function handles vector shape", {
  codegen_info <- list(dpar_shapes = list(shape = "vector"))
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    brms::brmsfamily("Gamma"), codegen_info
  )
  expect_match(code, "vector shape")
  expect_match(code, "shape\\[i\\], shape\\[i\\] / mu\\[i\\]")
})

test_that("subsampled beta function has correct call pattern", {
  codegen_info <- list(dpar_shapes = list(phi = "scalar"))
  code <- PDMPSamplersR:::make_subsampled_stan_functions(
    brms::brmsfamily("Beta"), codegen_info
  )
  expect_match(code, "subsampled_beta_lpdf")
  expect_match(code, "beta_lpdf")
  expect_match(code, "real phi")
  expect_match(code, "mu\\[i\\] \\* phi, \\(1 - mu\\[i\\]\\) \\* phi")
})

# -- Brace-aware block extraction ----------------------------------------------

test_that("extract_named_block extracts model block", {
  code <- paste(
    "model {",
    "  vector[N] mu = rep_vector(0.0, N);",
    "}",
    sep = "\n"
  )
  block <- PDMPSamplersR:::extract_named_block(code, "model")
  expect_match(block, "vector\\[N\\] mu")
})

test_that("extract_named_block is brace-aware with nested blocks", {
  code <- paste(
    "model {",
    "  if (1) {",
    "    vector[N] mu = rep_vector(0.0, N);",
    "  }",
    "}",
    "generated quantities {",
    "  for (n in 1:N) {",
    "    print(n);",
    "  }",
    "}",
    sep = "\n"
  )
  block <- PDMPSamplersR:::extract_named_block(code, "model")
  expect_match(block, "vector\\[N\\] mu")
  expect_no_match(block, "generated quantities")
})

test_that("extract_named_block errors on missing block", {
  code <- "data { int N; }"
  expect_error(
    PDMPSamplersR:::extract_named_block(code, "model"),
    "not found"
  )
})

test_that("replace_named_block replaces content correctly", {
  code <- "model {\n  old content\n}\ngenerated quantities {\n  gq\n}"
  result <- PDMPSamplersR:::replace_named_block(code, "model", "\n  new content\n")
  expect_match(result, "new content")
  expect_no_match(result, "old content")
  expect_match(result, "generated quantities")
})

# -- Rewrite helpers -----------------------------------------------------------

test_that("rewrite_mu_init resizes mu vector", {
  code <- "    vector[N] mu = rep_vector(0.0, N);"
  result <- PDMPSamplersR:::rewrite_mu_init(code)
  expect_match(result, "int m_sub = pdmp_get_subsample_size\\(\\)")
  expect_match(result, "vector\\[m_sub\\] mu = rep_vector\\(0\\.0, m_sub\\)")
  expect_no_match(result, "vector\\[N\\] mu")
})

test_that("rewrite_mu_init preserves indentation", {
  code <- "        vector[N] mu = rep_vector(0.0, N);"
  result <- PDMPSamplersR:::rewrite_mu_init(code)
  lines <- strsplit(result, "\n")[[1]]
  expect_true(all(grepl("^        ", lines)))
})

test_that("rewrite_mu_init errors on missing pattern", {
  code <- "  vector[N] eta = rep_vector(0.0, N);"
  expect_error(
    PDMPSamplersR:::rewrite_mu_init(code),
    "mu initialization"
  )
})

test_that("rewrite_fe_accumulation wraps Xc", {
  code <- "    mu += Intercept + Xc * b;"
  result <- PDMPSamplersR:::rewrite_fe_accumulation(code)
  expect_match(result, "get_subsampled_Xc\\(Xc\\) \\* b")
  expect_no_match(result, "(?<!subsampled_Xc\\()\\bXc\\b(?!\\))", perl = TRUE)
})

test_that("rewrite_fe_accumulation is no-op when Xc absent", {
  code <- "    mu += Intercept;"
  result <- PDMPSamplersR:::rewrite_fe_accumulation(code)
  expect_identical(result, code)
})

# -- Integration: validate_and_rewrite_subsampled_code -------------------------

test_that("validate_and_rewrite_subsampled_code rewrites a supported model", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    mu += Intercept + Xc * b;",
    "    target += subsampled_bernoulli_logit_lpmf(Y | mu, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(stancode)
  expect_match(result, "int m_sub = pdmp_get_subsample_size\\(\\)")
  expect_match(result, "vector\\[m_sub\\] mu")
  expect_match(result, "get_subsampled_Xc\\(Xc\\) \\* b")
  expect_no_match(result, "vector\\[N\\] mu")
})

# -- Integration: brm_stancode with new flow -----------------------------------

test_that("brm_stancode with subsampling returns custom-family code", {
  skip_if_no_brms()

  df <- data.frame(y = rbinom(50, 1, 0.5), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = brms::bernoulli(),
                          subsample_size = 10L)

  expect_true(is.list(result))
  expect_named(result, c("standard", "ext_cpp"))
  expect_match(result$ext_cpp, "subsampled_bernoulli_logit_lpmf")
  expect_match(result$ext_cpp, "pdmp_get_subsample_size")
  expect_match(result$ext_cpp, "get_subsampled_Xc")
  expect_match(result$ext_cpp, "int m_sub")
})

test_that("brm_stancode standard variant unchanged with new flow", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  original <- brms::stancode(y ~ x, data = df)
  result <- brm_stancode(y ~ x, data = df, subsample_size = 10L)

  expect_identical(result$standard, original)
})

test_that("brm_stancode with subsampling works for poisson", {
  skip_if_no_brms()

  df <- data.frame(y = rpois(50, 2), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = poisson(),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "subsampled_poisson_log_lpmf")
  expect_match(result$ext_cpp, "poisson_log_lpmf")
})

test_that("brm_stancode with subsampling works for gaussian", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = gaussian(),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "subsampled_gaussian_lpdf")
  expect_match(result$ext_cpp, "normal_lpdf")
  expect_match(result$ext_cpp, "real sigma")
})

test_that("brm_stancode accepts random effects with subsampling", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(50), x = rnorm(50), g = rep(1:5, each = 10))
  result <- brm_stancode(y ~ x + (1 | g), data = df, subsample_size = 10L)

  expect_true(is.list(result))
  expect_match(result$ext_cpp, "for \\(i in 1:m_sub\\)")
  expect_match(result$ext_cpp, "pdmp_get_subsample_index")
  expect_match(result$ext_cpp, "mu\\[i\\]")
})

test_that("brm_stancode accepts distributional formulas with subsampling", {
  skip_if_no_brms()

  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_stancode(brms::bf(y ~ x, sigma ~ x), data = df,
                          subsample_size = 10L)

  expect_true(is.list(result))
  expect_match(result$ext_cpp, "vector\\[m_sub\\] sigma")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xc_sigma\\)")
  expect_match(result$ext_cpp, "vector sigma")
  expect_match(result$ext_cpp, "sigma\\[i\\]")
})

# -- Snapshot tests for generated code -----------------------------------------

test_that("custom-family code snapshot: bernoulli fixed effects", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rbinom(50, 1, 0.5), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = brms::bernoulli(),
                          subsample_size = 10L)
  expect_snapshot(cat(result$ext_cpp))
})

test_that("custom-family code snapshot: poisson fixed effects", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rpois(50, 2), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = poisson(),
                          subsample_size = 10L)
  expect_snapshot(cat(result$ext_cpp))
})

test_that("custom-family code snapshot: gaussian fixed effects", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = gaussian(),
                          subsample_size = 10L)
  expect_snapshot(cat(result$ext_cpp))
})

# -- Spline rewrite tests -----------------------------------------------------

test_that("rewrite_spline_matrices wraps Xs", {
  code <- "    mu += Intercept + Xs * bs + Zs_1_1 * s_1_1;"
  result <- PDMPSamplersR:::rewrite_spline_matrices(code)
  expect_match(result, "get_subsampled_Xc\\(Xs\\) \\* bs")
  expect_match(result, "get_subsampled_Xc\\(Zs_1_1\\) \\* s_1_1")
})

test_that("rewrite_spline_matrices handles multiple spline bases", {
  code <- "    mu += Intercept + Xs * bs + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1;"
  result <- PDMPSamplersR:::rewrite_spline_matrices(code)
  expect_match(result, "get_subsampled_Xc\\(Zs_1_1\\)")
  expect_match(result, "get_subsampled_Xc\\(Zs_2_1\\)")
})

test_that("rewrite_spline_matrices is no-op without splines", {
  code <- "    mu += Intercept + Xc * b;"
  result <- PDMPSamplersR:::rewrite_spline_matrices(code)
  expect_identical(result, code)
})

# -- GP rewrite tests ----------------------------------------------------------

test_that("rewrite_gp_indexing wraps Jgp arrays", {
  code <- "    mu += Intercept + gp_pred_1[Jgp_1];"
  result <- PDMPSamplersR:::rewrite_gp_indexing(code)
  expect_match(result, "gp_pred_1\\[get_subsampled_int_array\\(Jgp_1\\)\\]")
})

test_that("rewrite_gp_indexing handles multiple GPs", {
  code <- "    mu += Intercept + gp_pred_1[Jgp_1] + gp_pred_2[Jgp_2];"
  result <- PDMPSamplersR:::rewrite_gp_indexing(code)
  expect_match(result, "gp_pred_1\\[get_subsampled_int_array\\(Jgp_1\\)\\]")
  expect_match(result, "gp_pred_2\\[get_subsampled_int_array\\(Jgp_2\\)\\]")
})

test_that("rewrite_gp_indexing is no-op without GPs", {
  code <- "    mu += Intercept + Xc * b;"
  result <- PDMPSamplersR:::rewrite_gp_indexing(code)
  expect_identical(result, code)
})

# -- Integration: spline model via validate_and_rewrite_subsampled_code --------

test_that("validate_and_rewrite handles spline-only model", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    mu += Intercept + Xs * bs + Zs_1_1 * s_1_1;",
    "    target += subsampled_gaussian_lpdf(Y | mu, sigma, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(stancode)
  expect_match(result, "int m_sub = pdmp_get_subsample_size\\(\\)")
  expect_match(result, "vector\\[m_sub\\] mu")
  expect_match(result, "get_subsampled_Xc\\(Xs\\) \\* bs")
  expect_match(result, "get_subsampled_Xc\\(Zs_1_1\\) \\* s_1_1")
})

test_that("validate_and_rewrite handles GP model", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[Nsubgp_1] gp_pred_1 = gp_exp_quad(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    mu += Intercept + gp_pred_1[Jgp_1];",
    "    target += subsampled_gaussian_lpdf(Y | mu, sigma, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(stancode)
  expect_match(result, "int m_sub = pdmp_get_subsample_size\\(\\)")
  expect_match(result, "gp_pred_1\\[get_subsampled_int_array\\(Jgp_1\\)\\]")
})

test_that("validate_and_rewrite handles FE + spline + GP combined", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[Nsubgp_1] gp_pred_1 = gp_exp_quad(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    mu += Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1 + gp_pred_1[Jgp_1];",
    "    target += subsampled_gaussian_lpdf(Y | mu, sigma, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(stancode)
  expect_match(result, "get_subsampled_Xc\\(Xc\\) \\* b")
  expect_match(result, "get_subsampled_Xc\\(Xs\\) \\* bs")
  expect_match(result, "get_subsampled_Xc\\(Zs_1_1\\) \\* s_1_1")
  expect_match(result, "gp_pred_1\\[get_subsampled_int_array\\(Jgp_1\\)\\]")
})

# -- Integration: brm_stancode with splines and GPs ---------------------------

test_that("brm_stancode with subsampling works for spline model", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_stancode(y ~ s(x), data = df, family = gaussian(),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xs\\)")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Zs_1_1\\)")
  expect_match(result$ext_cpp, "int m_sub")
})

test_that("brm_stancode with subsampling works for FE + spline model", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50), z = rnorm(50))
  result <- brm_stancode(y ~ x + s(z), data = df, family = gaussian(),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xc\\)")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xs\\)")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Zs_1_1\\)")
})

test_that("brm_stancode with subsampling works for GP model", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50))
  result <- brm_stancode(y ~ gp(x), data = df, family = gaussian(),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "get_subsampled_int_array\\(Jgp_1\\)")
  expect_match(result$ext_cpp, "int m_sub")
})

test_that("brm_stancode with subsampling works for bernoulli spline model", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rbinom(50, 1, 0.5), x = rnorm(50))
  result <- brm_stancode(y ~ s(x), data = df, family = brms::bernoulli(),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "subsampled_bernoulli_logit_lpmf")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xs\\)")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Zs_1_1\\)")
})

# -- Phase 2: Unified validation -----------------------------------------------

test_that("validate_subsampled_surface accepts FE only", {
  model_block <- paste(
    "  vector[N] mu = rep_vector(0.0, N);",
    "  mu += Intercept + Xc * b;",
    sep = "\n"
  )
  expect_no_error(PDMPSamplersR:::validate_subsampled_surface(model_block))
})

test_that("validate_subsampled_surface accepts RE + FE", {
  model_block <- paste(
    "  vector[N] mu = rep_vector(0.0, N);",
    "  mu += Intercept + Xc * b;",
    "  for (n in 1:N) {",
    "    mu[n] += r_1_1[J_1[n]] * Z_1_1[n];",
    "  }",
    sep = "\n"
  )
  expect_no_error(PDMPSamplersR:::validate_subsampled_surface(model_block))
})

test_that("validate_subsampled_surface accepts RE with multiple groups", {
  model_block <- paste(
    "  vector[N] mu = rep_vector(0.0, N);",
    "  mu += Intercept + Xc * b;",
    "  for (n in 1:N) {",
    "    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];",
    "  }",
    sep = "\n"
  )
  expect_no_error(PDMPSamplersR:::validate_subsampled_surface(model_block))
})

test_that("validate_subsampled_surface accepts RE with slope", {
  model_block <- paste(
    "  vector[N] mu = rep_vector(0.0, N);",
    "  mu += Intercept + Xc * b;",
    "  for (n in 1:N) {",
    "    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];",
    "  }",
    sep = "\n"
  )
  expect_no_error(PDMPSamplersR:::validate_subsampled_surface(model_block))
})

test_that("validate_subsampled_surface accepts FE + spline + RE", {
  model_block <- paste(
    "  vector[N] mu = rep_vector(0.0, N);",
    "  mu += Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1;",
    "  for (n in 1:N) {",
    "    mu[n] += r_1_1[J_1[n]] * Z_1_1[n];",
    "  }",
    sep = "\n"
  )
  expect_no_error(PDMPSamplersR:::validate_subsampled_surface(model_block))
})

test_that("validate_subsampled_surface rejects unknown dpar matrices", {
  model_block <- paste(
    "  vector[N] mu = rep_vector(0.0, N);",
    "  mu += Intercept + Xc * b;",
    "  vector[N] sigma = rep_vector(0.0, N);",
    "  sigma += Xc_sigma * b_sigma;",
    sep = "\n"
  )
  expect_error(
    PDMPSamplersR:::validate_subsampled_surface(model_block),
    "distributional predictor"
  )
})

test_that("validate_subsampled_surface accepts known dpar matrices", {
  model_block <- paste(
    "  vector[N] mu = rep_vector(0.0, N);",
    "  mu += Intercept + Xc * b;",
    "  vector[N] sigma = rep_vector(0.0, N);",
    "  sigma += Intercept_sigma + Xc_sigma * b_sigma;",
    sep = "\n"
  )
  expect_no_error(
    PDMPSamplersR:::validate_subsampled_surface(model_block, non_mu_dpars = "sigma")
  )
})

test_that("validate_re_loop_shape skips loop without mu[n]", {
  model_block <- paste(
    "  for (n in 1:N) {",
    "    eta[n] += r_1_1[J_1[n]] * Z_1_1[n];",
    "  }",
    sep = "\n"
  )
  expect_no_error(PDMPSamplersR:::validate_re_loop_shape(model_block))
})

test_that("validate_re_loop_shape rejects loop without r_G_K pattern", {
  model_block <- paste(
    "  for (n in 1:N) {",
    "    mu[n] += some_other_thing[n];",
    "  }",
    sep = "\n"
  )
  expect_error(
    PDMPSamplersR:::validate_re_loop_shape(model_block),
    "r_G_K"
  )
})

# -- Phase 2: RE loop rewriting ------------------------------------------------

test_that("rewrite_re_loops rewrites single-group intercept loop", {
  code <- paste(
    "    for (n in 1:N) {",
    "      // add more terms to the linear predictor",
    "      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];",
    "    }",
    sep = "\n"
  )
  result <- PDMPSamplersR:::rewrite_re_loops(code)
  expect_match(result, "for \\(i in 1:m_sub\\)")
  expect_match(result, "int n = pdmp_get_subsample_index\\(i\\);")
  expect_match(result, "mu\\[i\\]")
  expect_no_match(result, "mu\\[n\\]")
  expect_match(result, "r_1_1\\[J_1\\[n\\]\\]")
  expect_match(result, "Z_1_1\\[n\\]")
})

test_that("rewrite_re_loops rewrites multi-group loop", {
  code <- paste(
    "    for (n in 1:N) {",
    "      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];",
    "    }",
    sep = "\n"
  )
  result <- PDMPSamplersR:::rewrite_re_loops(code)
  expect_match(result, "for \\(i in 1:m_sub\\)")
  expect_match(result, "mu\\[i\\]")
  expect_match(result, "r_1_1\\[J_1\\[n\\]\\]")
  expect_match(result, "r_2_1\\[J_2\\[n\\]\\]")
})

test_that("rewrite_re_loops rewrites slope loop", {
  code <- paste(
    "    for (n in 1:N) {",
    "      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];",
    "    }",
    sep = "\n"
  )
  result <- PDMPSamplersR:::rewrite_re_loops(code)
  expect_match(result, "for \\(i in 1:m_sub\\)")
  expect_match(result, "mu\\[i\\]")
  expect_match(result, "r_1_2\\[J_1\\[n\\]\\]")
  expect_match(result, "Z_1_2\\[n\\]")
})

test_that("rewrite_re_loops is no-op without RE loop", {
  code <- "    mu += Intercept + Xc * b;"
  result <- PDMPSamplersR:::rewrite_re_loops(code)
  expect_identical(result, code)
})

# -- Phase 2 integration: validate_and_rewrite with RE ------------------------

test_that("validate_and_rewrite handles FE + RE model", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    mu += Intercept + Xc * b;",
    "    for (n in 1:N) {",
    "      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];",
    "    }",
    "    target += subsampled_bernoulli_logit_lpmf(Y | mu, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(stancode)
  expect_match(result, "int m_sub = pdmp_get_subsample_size\\(\\)")
  expect_match(result, "vector\\[m_sub\\] mu")
  expect_match(result, "get_subsampled_Xc\\(Xc\\) \\* b")
  expect_match(result, "for \\(i in 1:m_sub\\)")
  expect_match(result, "int n = pdmp_get_subsample_index\\(i\\);")
  expect_match(result, "mu\\[i\\]")
  expect_no_match(result, "for \\(n in 1:N\\)")
  expect_no_match(result, "mu\\[n\\]")
})

test_that("validate_and_rewrite handles FE + spline + RE model", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    mu += Intercept + Xc * b + Xs * bs + Zs_1_1 * s_1_1;",
    "    for (n in 1:N) {",
    "      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];",
    "    }",
    "    target += subsampled_gaussian_lpdf(Y | mu, sigma, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(stancode)
  expect_match(result, "get_subsampled_Xc\\(Xc\\) \\* b")
  expect_match(result, "get_subsampled_Xc\\(Xs\\) \\* bs")
  expect_match(result, "get_subsampled_Xc\\(Zs_1_1\\) \\* s_1_1")
  expect_match(result, "for \\(i in 1:m_sub\\)")
  expect_match(result, "mu\\[i\\]")
})

# -- Phase 2 integration: brm_stancode with RE --------------------------------

test_that("brm_stancode with subsampling works for RE intercept model", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50), g = rep(1:5, each = 10))
  result <- brm_stancode(y ~ x + (1 | g), data = df, family = gaussian(),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "for \\(i in 1:m_sub\\)")
  expect_match(result$ext_cpp, "int n = pdmp_get_subsample_index\\(i\\);")
  expect_match(result$ext_cpp, "mu\\[i\\]")
  expect_match(result$ext_cpp, "r_1_1\\[J_1\\[n\\]\\]")
})

test_that("brm_stancode with subsampling works for RE slope model", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50), g = rep(1:5, each = 10))
  result <- brm_stancode(y ~ x + (1 + x | g), data = df, family = gaussian(),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "for \\(i in 1:m_sub\\)")
  expect_match(result$ext_cpp, "r_1_1\\[J_1\\[n\\]\\]")
  expect_match(result$ext_cpp, "r_1_2\\[J_1\\[n\\]\\]")
})

test_that("brm_stancode with subsampling works for multi-group RE model", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rbinom(50, 1, 0.5), x = rnorm(50),
                   g1 = rep(1:5, each = 10), g2 = rep(1:10, each = 5))
  result <- brm_stancode(y ~ x + (1 | g1) + (1 | g2), data = df,
                          family = brms::bernoulli(), subsample_size = 10L)

  expect_match(result$ext_cpp, "for \\(i in 1:m_sub\\)")
  expect_match(result$ext_cpp, "r_1_1\\[J_1\\[n\\]\\]")
  expect_match(result$ext_cpp, "r_2_1\\[J_2\\[n\\]\\]")
})

test_that("brm_stancode with subsampling works for FE + spline + RE", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(60), x = rnorm(60),
                   z = seq(0, 2 * pi, length.out = 60),
                   g = rep(1:6, each = 10))
  result <- brm_stancode(y ~ x + s(z, k = 5) + (1 | g), data = df,
                          family = gaussian(), subsample_size = 10L)

  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xc\\)")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xs\\)")
  expect_match(result$ext_cpp, "for \\(i in 1:m_sub\\)")
  expect_match(result$ext_cpp, "r_1_1\\[J_1\\[n\\]\\]")
})

# -- Phase 3: Dpar rewrite unit tests ------------------------------------------

test_that("rewrite_dpar_init resizes dpar vector", {
  code <- "    vector[N] sigma = rep_vector(0.0, N);"
  result <- PDMPSamplersR:::rewrite_dpar_init(code, "sigma")
  expect_match(result, "vector\\[m_sub\\] sigma = rep_vector\\(0\\.0, m_sub\\)")
  expect_no_match(result, "vector\\[N\\] sigma")
})

test_that("rewrite_dpar_init works for shape", {
  code <- "    vector[N] shape = rep_vector(0.0, N);"
  result <- PDMPSamplersR:::rewrite_dpar_init(code, "shape")
  expect_match(result, "vector\\[m_sub\\] shape = rep_vector\\(0\\.0, m_sub\\)")
})

test_that("rewrite_dpar_fe wraps Xc_sigma", {
  code <- "    sigma += Intercept_sigma + Xc_sigma * b_sigma;"
  result <- PDMPSamplersR:::rewrite_dpar_fe(code, "sigma")
  expect_match(result, "get_subsampled_Xc\\(Xc_sigma\\) \\* b_sigma")
})

test_that("rewrite_dpar_fe wraps Xc_shape", {
  code <- "    shape += Intercept_shape + Xc_shape * b_shape;"
  result <- PDMPSamplersR:::rewrite_dpar_fe(code, "shape")
  expect_match(result, "get_subsampled_Xc\\(Xc_shape\\) \\* b_shape")
})

test_that("rewrite_dpar_fe is no-op without dpar FE", {
  code <- "    sigma += Intercept_sigma;"
  result <- PDMPSamplersR:::rewrite_dpar_fe(code, "sigma")
  expect_identical(result, code)
})

test_that("rewrite_dpar_splines wraps Xs_sigma and Zs_sigma_*", {
  code <- "    sigma += Intercept_sigma + Xs_sigma * bs_sigma + Zs_sigma_1_1 * s_sigma_1_1;"
  result <- PDMPSamplersR:::rewrite_dpar_splines(code, "sigma")
  expect_match(result, "get_subsampled_Xc\\(Xs_sigma\\) \\* bs_sigma")
  expect_match(result, "get_subsampled_Xc\\(Zs_sigma_1_1\\) \\* s_sigma_1_1")
})

test_that("rewrite_dpar_splines is no-op without dpar splines", {
  code <- "    sigma += Intercept_sigma + Xc_sigma * b_sigma;"
  result <- PDMPSamplersR:::rewrite_dpar_splines(code, "sigma")
  expect_identical(result, code)
})

test_that("rewrite_dpar_gp_indexing wraps Jgp_sigma arrays", {
  code <- "    sigma += Intercept_sigma + gp_pred_sigma_1[Jgp_sigma_1];"
  result <- PDMPSamplersR:::rewrite_dpar_gp_indexing(code, "sigma")
  expect_match(result, "gp_pred_sigma_1\\[get_subsampled_int_array\\(Jgp_sigma_1\\)\\]")
})

test_that("rewrite_dpar_gp_indexing is no-op without dpar GPs", {
  code <- "    sigma += Intercept_sigma + Xc_sigma * b_sigma;"
  result <- PDMPSamplersR:::rewrite_dpar_gp_indexing(code, "sigma")
  expect_identical(result, code)
})

test_that("rewrite_re_loops rewrites dpar accumulations", {
  code <- paste(
    "    for (n in 1:N) {",
    "      sigma[n] += r_1_sigma_1[J_1[n]] * Z_1_sigma_1[n];",
    "    }",
    sep = "\n"
  )
  result <- PDMPSamplersR:::rewrite_re_loops(code, non_mu_dpars = "sigma")
  expect_match(result, "for \\(i in 1:m_sub\\)")
  expect_match(result, "sigma\\[i\\]")
  expect_no_match(result, "sigma\\[n\\]")
})

test_that("rewrite_re_loops rewrites both mu and dpar loops", {
  code <- paste(
    "    for (n in 1:N) {",
    "      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];",
    "    }",
    "    for (n in 1:N) {",
    "      sigma[n] += r_2_sigma_1[J_2[n]] * Z_2_sigma_1[n];",
    "    }",
    sep = "\n"
  )
  result <- PDMPSamplersR:::rewrite_re_loops(code, non_mu_dpars = "sigma")
  expect_match(result, "mu\\[i\\]")
  expect_match(result, "sigma\\[i\\]")
  expect_no_match(result, "mu\\[n\\]")
  expect_no_match(result, "sigma\\[n\\]")
})

# -- Phase 3: validate_subsampled_surface with dpar RE -------------------------

test_that("validate_subsampled_surface accepts dpar RE loops", {
  model_block <- paste(
    "  vector[N] mu = rep_vector(0.0, N);",
    "  mu += Intercept + Xc * b;",
    "  vector[N] sigma = rep_vector(0.0, N);",
    "  sigma += Intercept_sigma;",
    "  for (n in 1:N) {",
    "    sigma[n] += r_1_sigma_1[J_1[n]] * Z_1_sigma_1[n];",
    "  }",
    sep = "\n"
  )
  expect_no_error(
    PDMPSamplersR:::validate_subsampled_surface(model_block, non_mu_dpars = "sigma")
  )
})

test_that("validate_subsampled_surface accepts combined mu+dpar RE", {
  model_block <- paste(
    "  vector[N] mu = rep_vector(0.0, N);",
    "  mu += Intercept + Xc * b;",
    "  vector[N] sigma = rep_vector(0.0, N);",
    "  sigma += Intercept_sigma;",
    "  for (n in 1:N) {",
    "    mu[n] += r_1_1[J_1[n]] * Z_1_1[n];",
    "  }",
    "  for (n in 1:N) {",
    "    sigma[n] += r_2_sigma_1[J_2[n]] * Z_2_sigma_1[n];",
    "  }",
    sep = "\n"
  )
  expect_no_error(
    PDMPSamplersR:::validate_subsampled_surface(model_block, non_mu_dpars = "sigma")
  )
})

# -- Phase 3 integration: validate_and_rewrite with dpar ----------------------

test_that("validate_and_rewrite handles dpar FE model", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    vector[N] sigma = rep_vector(0.0, N);",
    "    mu += Intercept + Xc * b;",
    "    sigma += Intercept_sigma + Xc_sigma * b_sigma;",
    "    sigma = exp(sigma);",
    "    target += subsampled_gaussian_lpdf(Y | mu, sigma, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(
    stancode, non_mu_dpars = "sigma"
  )
  expect_match(result, "vector\\[m_sub\\] mu")
  expect_match(result, "vector\\[m_sub\\] sigma")
  expect_match(result, "get_subsampled_Xc\\(Xc\\) \\* b")
  expect_match(result, "get_subsampled_Xc\\(Xc_sigma\\) \\* b_sigma")
  expect_no_match(result, "vector\\[N\\] sigma")
})

test_that("validate_and_rewrite handles dpar spline model", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    vector[N] sigma = rep_vector(0.0, N);",
    "    mu += Intercept + Xc * b;",
    "    sigma += Intercept_sigma + Xs_sigma * bs_sigma + Zs_sigma_1_1 * s_sigma_1_1;",
    "    sigma = exp(sigma);",
    "    target += subsampled_gaussian_lpdf(Y | mu, sigma, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(
    stancode, non_mu_dpars = "sigma"
  )
  expect_match(result, "vector\\[m_sub\\] sigma")
  expect_match(result, "get_subsampled_Xc\\(Xs_sigma\\) \\* bs_sigma")
  expect_match(result, "get_subsampled_Xc\\(Zs_sigma_1_1\\) \\* s_sigma_1_1")
})

test_that("validate_and_rewrite handles dpar GP model", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[Nsubgp_sigma_1] gp_pred_sigma_1 = gp_exp_quad(Xgp_sigma_1, sdgp_sigma_1[1], lscale_sigma_1[1], zgp_sigma_1);",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    vector[N] sigma = rep_vector(0.0, N);",
    "    mu += Intercept + Xc * b;",
    "    sigma += Intercept_sigma + gp_pred_sigma_1[Jgp_sigma_1];",
    "    sigma = exp(sigma);",
    "    target += subsampled_gaussian_lpdf(Y | mu, sigma, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(
    stancode, non_mu_dpars = "sigma"
  )
  expect_match(result, "vector\\[m_sub\\] sigma")
  expect_match(result, "gp_pred_sigma_1\\[get_subsampled_int_array\\(Jgp_sigma_1\\)\\]")
})

test_that("validate_and_rewrite handles dpar RE model", {
  stancode <- paste(
    "functions {",
    "  int pdmp_get_subsample_size();",
    "}",
    "data {",
    "  int<lower=1> N;",
    "}",
    "model {",
    "  if (!prior_only) {",
    "    vector[N] mu = rep_vector(0.0, N);",
    "    vector[N] sigma = rep_vector(0.0, N);",
    "    mu += Intercept + Xc * b;",
    "    sigma += Intercept_sigma;",
    "    for (n in 1:N) {",
    "      sigma[n] += r_1_sigma_1[J_1[n]] * Z_1_sigma_1[n];",
    "    }",
    "    sigma = exp(sigma);",
    "    target += subsampled_gaussian_lpdf(Y | mu, sigma, N);",
    "  }",
    "}",
    sep = "\n"
  )
  result <- PDMPSamplersR:::validate_and_rewrite_subsampled_code(
    stancode, non_mu_dpars = "sigma"
  )
  expect_match(result, "vector\\[m_sub\\] sigma")
  expect_match(result, "for \\(i in 1:m_sub\\)")
  expect_match(result, "sigma\\[i\\]")
  expect_no_match(result, "sigma\\[n\\]")
})

# -- Phase 3 integration: brm_stancode with distributional formulas ------------

test_that("brm_stancode works for gaussian bf(y ~ x, sigma ~ z)", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50), z = rnorm(50))
  result <- brm_stancode(brms::bf(y ~ x, sigma ~ z), data = df,
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "vector\\[m_sub\\] sigma")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xc_sigma\\)")
  expect_match(result$ext_cpp, "vector sigma")
  expect_match(result$ext_cpp, "sigma\\[i\\]")
  expect_match(result$ext_cpp, "sigma = exp\\(sigma\\)")
})

test_that("brm_stancode works for negbinomial bf(y ~ x, shape ~ z)", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnbinom(50, size = 5, mu = 2), x = rnorm(50), z = rnorm(50))
  result <- brm_stancode(brms::bf(y ~ x, shape ~ z), data = df,
                          family = brms::negbinomial(), subsample_size = 10L)

  expect_match(result$ext_cpp, "vector\\[m_sub\\] shape")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xc_shape\\)")
  expect_match(result$ext_cpp, "vector shape")
  expect_match(result$ext_cpp, "shape\\[i\\]")
})

test_that("brm_stancode works for gaussian sigma ~ s(z)", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(100), x = rnorm(100),
                   z = seq(0, 2 * pi, length.out = 100))
  result <- brm_stancode(brms::bf(y ~ x, sigma ~ s(z, k = 5)), data = df,
                          subsample_size = 20L)

  expect_match(result$ext_cpp, "vector\\[m_sub\\] sigma")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xs_sigma\\)")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Zs_sigma_1_1\\)")
})

test_that("brm_stancode works for gaussian sigma ~ 1 + (1|g)", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rnorm(50), x = rnorm(50), g = rep(1:5, each = 10))
  result <- brm_stancode(brms::bf(y ~ x, sigma ~ 1 + (1 | g)), data = df,
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "vector\\[m_sub\\] sigma")
  expect_match(result$ext_cpp, "for \\(i in 1:m_sub\\)")
  expect_match(result$ext_cpp, "sigma\\[i\\]")
})

# -- Phase 4: brm_stancode with new families -----------------------------------

test_that("brm_stancode works for student y ~ x", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = brms::rstudent_t(50, df = 3, mu = 1, sigma = 2),
                   x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = brms::student(),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "subsampled_student_lpdf")
  expect_match(result$ext_cpp, "student_t_lpdf")
  expect_match(result$ext_cpp, "nu, mu\\[i\\], sigma")
  expect_match(result$ext_cpp, "logm1")
  expect_match(result$ext_cpp, "expp1")
})

test_that("brm_stancode works for Gamma y ~ x", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rgamma(50, shape = 2, rate = 1), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = brms::brmsfamily("Gamma"),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "subsampled_gamma_lpdf")
  expect_match(result$ext_cpp, "gamma_lpdf")
  expect_match(result$ext_cpp, "shape, shape / mu\\[i\\]")
  expect_match(result$ext_cpp, "mu = exp\\(mu\\)")
})

test_that("brm_stancode works for Beta y ~ x", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rbeta(50, 2, 3), x = rnorm(50))
  result <- brm_stancode(y ~ x, data = df, family = brms::brmsfamily("Beta"),
                          subsample_size = 10L)

  expect_match(result$ext_cpp, "subsampled_beta_lpdf")
  expect_match(result$ext_cpp, "beta_lpdf")
  expect_match(result$ext_cpp, "mu\\[i\\] \\* phi")
  expect_match(result$ext_cpp, "mu = inv_logit\\(mu\\)")
})

test_that("brm_stancode works for student bf(y ~ x, sigma ~ z)", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = brms::rstudent_t(50, df = 3, mu = 1, sigma = 2),
                   x = rnorm(50), z = rnorm(50))
  result <- brm_stancode(brms::bf(y ~ x, sigma ~ z), data = df,
                          family = brms::student(), subsample_size = 10L)

  expect_match(result$ext_cpp, "vector\\[m_sub\\] sigma")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xc_sigma\\)")
  expect_match(result$ext_cpp, "sigma\\[i\\]")
})

test_that("brm_stancode works for Gamma bf(y ~ x, shape ~ z)", {
  skip_if_no_brms()

  set.seed(42)
  df <- data.frame(y = rgamma(50, shape = 2, rate = 1),
                   x = rnorm(50), z = rnorm(50))
  result <- brm_stancode(brms::bf(y ~ x, shape ~ z), data = df,
                          family = brms::brmsfamily("Gamma"), subsample_size = 10L)

  expect_match(result$ext_cpp, "vector\\[m_sub\\] shape")
  expect_match(result$ext_cpp, "get_subsampled_Xc\\(Xc_shape\\)")
  expect_match(result$ext_cpp, "shape\\[i\\]")
})
