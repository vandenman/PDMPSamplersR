# Extracted from test-brm_stancode_standata.R:164

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
sdata <- brms::standata(y ~ x, data = df)
expected_means <- colMeans(sdata$X[, -1, drop = FALSE])
expect_equal(as.numeric(result$full$means_X), as.numeric(expected_means))
expect_equal(as.numeric(result$prior$means_X), as.numeric(expected_means))
expect_equal(as.numeric(result$subsample$means_X), as.numeric(expected_means))
