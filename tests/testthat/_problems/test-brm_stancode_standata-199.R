# Extracted from test-brm_stancode_standata.R:199

# prequel ----------------------------------------------------------------------
skip_if_no_brms <- function() {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("rstan")
}

# test -------------------------------------------------------------------------
skip_if_no_brms()
df <- data.frame(y = rnorm(50), x = rnorm(50))
code <- brm_stancode(y ~ x, data = df, subsample_size = 10L)
dlist <- brm_standata(y ~ x, data = df, subsample_size = 10L)
expect_true(grepl("means_X", code$standard))
expect_true("means_X" %in% names(dlist$full))
