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

# run these immediately so they are cached
is_julia_available()
is_pdmp_julia_backend_available()

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

brms_backend <- function() {
  backend <- Sys.getenv("PDMPSAMPLERSR_BRMS_BACKEND", "rstan")
  if (!backend %in% c("rstan", "cmdstanr")) {
    stop("Unsupported PDMPSAMPLERSR_BRMS_BACKEND: ", backend, call. = FALSE)
  }
  backend
}

skip_if_no_brms_backend <- function() {
  if (identical(brms_backend(), "cmdstanr")) {
    testthat::skip_if_not_installed("cmdstanr")

    cmdstan_available <- cache_get_or_set("cmdstan_available", tryCatch({
      !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))
    }, error = function(e) FALSE))

    testthat::skip_if_not(cmdstan_available, "CmdStan is not available")
  } else {
    testthat::skip_if_not_installed("rstan")
  }
}

skip_if_no_brms_setup <- function() {
  testthat::skip_if_not_installed("brms")
  skip_if_no_brms_backend()

  bridge_available <- cache_get_or_set("brms_bridge_setup_available", tryCatch({
    PDMPSamplersR:::check_for_julia_setup()
    # Check if the BridgeStan extension was triggered and thus this method is available.
    isTRUE(JuliaCall::julia_eval("hasmethod(PDMPModel, Tuple{String, String})"))
  }, error = function(e) FALSE))

  testthat::skip_if_not(bridge_available, "Julia + BridgeStan not available")
}
