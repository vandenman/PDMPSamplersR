# Extracted from test-brm_stancode_standata.R:85

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
expect_equal(result$full$N, 50L)
expect_true("means_X" %in% names(result$full))
expect_equal(length(result$full$means_X), result$full$Kc)
