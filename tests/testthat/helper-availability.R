.test_availability_cache <- new.env(parent = emptyenv())

cache_get_or_set <- function(key, expr) {
  if (!exists(key, envir = .test_availability_cache, inherits = FALSE)) {
    assign(key, expr, envir = .test_availability_cache)
  }
  get(key, envir = .test_availability_cache, inherits = FALSE)
}

is_julia_available <- function() {
  cache_get_or_set("julia_available", tryCatch({
    JuliaCall::julia_setup(verbose = FALSE)
    TRUE
  }, error = function(e) FALSE))
}

is_pdmp_julia_backend_available <- function() {
  cache_get_or_set("pdmp_julia_backend_available", tryCatch({
    PDMPSamplersR:::check_for_julia_setup()
    TRUE
  }, error = function(e) FALSE))
}

skip_if_no_julia <- function() {
  testthat::skip_if_not(is_julia_available(), "Julia is not available")
}

skip_if_no_pdmp_julia_backend <- function() {
  skip_if_no_julia()
  testthat::skip_if_not(
    is_pdmp_julia_backend_available(),
    "PDMPSamplers.jl Julia integration not functional"
  )
}

skip_if_no_brms_setup <- function() {
  testthat::skip_if_not_installed("brms")
  testthat::skip_if_not_installed("rstan")

  bridge_available <- cache_get_or_set("brms_bridge_setup_available", tryCatch({
    PDMPSamplersR:::check_for_julia_setup()
    isTRUE(JuliaCall::julia_eval("hasmethod(PDMPModel, Tuple{String, String})"))
  }, error = function(e) FALSE))

  testthat::skip_if_not(bridge_available, "Julia + BridgeStan not available")
}