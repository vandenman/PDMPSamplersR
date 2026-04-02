# Extracted from test-brm_stancode_standata.R:97

# prequel ----------------------------------------------------------------------
skip_if_no_brms <- function() {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("rstan")
}

# test -------------------------------------------------------------------------
skip_if_no_brms()
set.seed(42)
df <- data.frame(y = rnorm(50), x = rnorm(50))
result <- brm_standata(y ~ x, data = df, subsample_size = 10L)
expect_equal(result$prior$N, 1L)
expect_equal(result$prior$prior_only, 1L)
expect_true("means_X" %in% names(result$prior))
