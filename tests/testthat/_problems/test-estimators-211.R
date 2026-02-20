# Extracted from test-estimators.R:211

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "PDMPSamplersR", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
skip_if_no_julia <- function() {
  julia_available <- tryCatch({
    JuliaCall::julia_setup(verbose = FALSE)
    TRUE
  }, error = function(e) FALSE)
  testthat::skip_if_not(julia_available, "Julia is not available")
}

# test -------------------------------------------------------------------------
skip_on_cran()
skip_if_no_julia()
d <- 2
neg_grad <- function(x) x
result <- pdmp_sample(neg_grad, d = d, flow = "ZigZag", T = 2000,
                        show_progress = FALSE)
transforms <- list(identity_transform(), identity_transform())
m_identity <- mean(result, transforms = transforms)
