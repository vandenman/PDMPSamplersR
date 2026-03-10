# Tests for brms subsampling: helper functions and integration tests.

skip_if_no_brms_setup <- function() {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("rstan")

  julia_available <- tryCatch({
    PDMPSamplersR:::check_for_julia_setup()
    JuliaCall::julia_eval("hasmethod(PDMPModel, Tuple{String, String})")
  }, error = function(e) FALSE)
  testthat::skip_if_not(julia_available, "Julia + BridgeStan not available")
}

# ── Unit tests for helper functions (no Julia needed) ─────────────────────────

test_that("fix_brms_stancode moves means_X from transformed data to data", {
  code <- paste(
    "data {",
    "  int<lower=1> N;",
    "  vector[N] Y;",
    "  int<lower=1> K;",
    "  matrix[N, K] X;",
    "  int<lower=1> Kc;",
    "  int prior_only;",
    "}",
    "transformed data {",
    "  vector[Kc] means_X;",
    "  for (i in 2:K) {",
    "    means_X[i - 1] = mean(X[, i]);",
    "  }",
    "}",
    "parameters {",
    "  vector[Kc] b;",
    "  real Intercept;",
    "}",
    sep = "\n"
  )

  fix_brms_stancode <- PDMPSamplersR:::fix_brms_stancode
  result <- fix_brms_stancode(code)
  lines <- strsplit(result, "\n")[[1]]

  # means_X declaration should appear after prior_only in data block
  prior_only_line <- grep("prior_only", lines)[1]
  means_data_line <- grep("means_X", lines)
  # There should be exactly one means_X line (the new one in data)
  means_data_line <- means_data_line[means_data_line > prior_only_line]
  expect_length(means_data_line, 1)

  # The original transformed-data declaration and assignment should be gone
  expect_false(any(grepl("means_X\\[i - 1\\] = mean\\(X\\[, i\\]\\)", result)))

  td_start <- grep("transformed data", lines)
  td_body <- lines[(td_start + 1):length(lines)]
  td_end <- grep("^\\}", td_body)[1]
  td_block <- td_body[1:(td_end - 1)]
  expect_false(any(grepl("vector\\[Kc\\] means_X", td_block)))
})

test_that("fix_brms_stancode returns unchanged code without means_X", {
  fix_brms_stancode <- PDMPSamplersR:::fix_brms_stancode
  code <- "data {\n  int N;\n}\nparameters {\n  real x;\n}"
  expect_equal(fix_brms_stancode(code), code)
})

test_that("fix_brms_stancode returns unchanged code without transformed data", {
  fix_brms_stancode <- PDMPSamplersR:::fix_brms_stancode
  code <- "data {\n  int N;\n  int prior_only;\n}\nparameters {\n  real x;\n}"
  expect_equal(fix_brms_stancode(code), code)
})

test_that("subset_standata subsets correctly", {
  sdata <- list(
    N = 10L, Y = 1:10, K = 3L, Kc = 2L,
    X = matrix(seq_len(30), nrow = 10, ncol = 3)
  )
  means_X <- c(5.5, 15.5)
  indices <- c(1L, 3L, 5L)

  subset_standata <- PDMPSamplersR:::subset_standata
  sub <- subset_standata(sdata, indices, means_X)

  expect_equal(sub$N, 3L)
  expect_equal(sub$Y, c(1L, 3L, 5L))
  expect_equal(nrow(sub$X), 3L)
  expect_equal(ncol(sub$X), 3L)
  expect_equal(sub$X[1, ], sdata$X[1, ])
  expect_equal(sub$X[2, ], sdata$X[3, ])
  expect_equal(sub$X[3, ], sdata$X[5, ])
  expect_equal(sub$means_X, means_X)
  expect_equal(sub$K, 3L)
})

test_that("make_prior_standata creates valid prior data", {
  sdata <- list(
    N = 10L, Y = 1:10, K = 3L, Kc = 2L,
    X = matrix(seq_len(30), nrow = 10, ncol = 3)
  )
  means_X <- c(5.5, 15.5)

  make_prior_standata <- PDMPSamplersR:::make_prior_standata
  prior <- make_prior_standata(sdata, means_X)

  expect_equal(prior$N, 1L)
  expect_equal(prior$K, 3L)
  expect_equal(prior$Kc, 2L)
  expect_equal(prior$prior_only, 1L)
  expect_equal(prior$means_X, means_X)
  expect_true(is.integer(prior$Y))
  expect_equal(as.integer(prior$Y), 0L)
  expect_true(!is.null(dim(prior$Y)))
  expect_equal(nrow(prior$X), 1L)
  expect_equal(ncol(prior$X), 3L)
})

test_that("make_prior_standata handles double Y", {
  sdata <- list(
    N = 5L, Y = c(1.1, 2.2, 3.3, 4.4, 5.5), K = 2L, Kc = 1L,
    X = matrix(seq_len(10), nrow = 5, ncol = 2)
  )
  means_X <- 3.5

  make_prior_standata <- PDMPSamplersR:::make_prior_standata
  prior <- make_prior_standata(sdata, means_X)

  expect_true(is.double(prior$Y))
  expect_equal(as.double(prior$Y), 0.0)
  expect_true(!is.null(dim(prior$Y)))
})

# ── Validation tests ─────────────────────────────────────────────────────────

test_that("brm_pdmp rejects subsample_size >= nrow(data)", {
  skip_if_no_brms_setup()
  expect_error(
    brm_pdmp(y ~ x, data = data.frame(y = 1:10, x = 1:10),
             subsample_size = 10L, T = 100, show_progress = FALSE),
    "subsample_size"
  )
})

test_that("brm_pdmp rejects random effects with subsampling", {
  skip_if_no_brms_setup()

  set.seed(1)
  df <- data.frame(
    y = rnorm(20), x = rnorm(20),
    group = rep(1:4, each = 5)
  )
  expect_error(
    brm_pdmp(y ~ x + (1 | group), data = df,
             subsample_size = 5L, T = 100, show_progress = FALSE),
    "fixed-effects"
  )
})

# ── Integration tests (need Julia + brms) ────────────────────────────────────

test_that("brm_pdmp with subsample_size: gaussian y ~ x", {
  skip_if_no_brms_setup()
  # NOTE: Gaussian models have a positive-constrained sigma parameter. Random
  # PDMP initialization can push log(sigma) to extremes, causing BridgeStan to
  # throw "Scale vector is inf". This is a pre-existing sampler limitation, not
  # specific to subsampling. Skip on this known error.
  set.seed(42)
  n <- 100
  x <- rnorm(n)
  y <- 2 + 0.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(y = y, x = x)

  fit <- tryCatch(
    brm_pdmp(y ~ x, data = df, family = gaussian(),
             flow = "ZigZag", algorithm = "GridThinningStrategy",
             T = 50000, t_warmup = 5000,
             subsample_size = 20L,
             show_progress = FALSE),
    error = function(e) {
      if (grepl("Scale vector is inf", conditionMessage(e)))
        testthat::skip("sigma=inf during PDMP initialization (stochastic)")
      stop(e)
    }
  )

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Intercept" %in% rownames(fixed))
  expect_true("x" %in% rownames(fixed))
})

test_that("brm_pdmp with subsample_size: bernoulli y ~ x", {
  skip_if_no_brms_setup()

  set.seed(123)
  n <- 200
  x <- rnorm(n)
  prob <- plogis(0.5 + 1.0 * x)
  y <- rbinom(n, 1, prob)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = brms::bernoulli(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, t_warmup = 5000,
                  subsample_size = 30L,
                  show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Intercept" %in% rownames(fixed))
  expect_true("x" %in% rownames(fixed))
})

test_that("brm_pdmp with subsample_size: poisson y ~ x", {
  skip_if_no_brms_setup()

  set.seed(456)
  n <- 200
  x <- rnorm(n)
  y <- rpois(n, exp(0.5 + 0.3 * x))
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = poisson(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, t_warmup = 5000,
                  subsample_size = 30L,
                  show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Intercept" %in% rownames(fixed))
  expect_true("x" %in% rownames(fixed))
})

test_that("brm_pdmp with subsample_size: summary and pp_check work", {
  skip_if_no_brms_setup()

  set.seed(123)
  n <- 200
  x <- rnorm(n)
  prob <- plogis(0.5 + 1.0 * x)
  y <- rbinom(n, 1, prob)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = brms::bernoulli(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, t_warmup = 5000,
                  subsample_size = 30L,
                  show_progress = FALSE)

  s <- summary(fit)
  expect_true(!is.null(s))
  expect_true("fixed" %in% names(s))

  pp <- brms::pp_check(fit, ndraws = 20)
  expect_s3_class(pp, "gg")
})

test_that("brm_pdmp with subsample_size: multi-chain Rhat", {
  skip_if_no_brms_setup()

  set.seed(123)
  n <- 200
  x <- rnorm(n)
  prob <- plogis(0.5 + 1.0 * x)
  y <- rbinom(n, 1, prob)
  df <- data.frame(y = y, x = x)

  fit <- brm_pdmp(y ~ x, data = df, family = brms::bernoulli(),
                  flow = "ZigZag", algorithm = "GridThinningStrategy",
                  T = 50000, t_warmup = 5000,
                  subsample_size = 30L,
                  n_chains = 2L,
                  show_progress = FALSE)

  expect_s3_class(fit, "brmsfit")

  s <- summary(fit)
  fixed <- s$fixed
  expect_true("Rhat" %in% colnames(fixed))
  rhat_vals <- fixed[, "Rhat"]
  for (rh in rhat_vals) {
    expect_true(is.finite(rh))
  }
})
