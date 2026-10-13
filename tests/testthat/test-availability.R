test_that("the no-Julia job never initializes Julia", {
  old <- Sys.getenv("PDMPSAMPLERSR_TEST_JULIA", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("PDMPSAMPLERSR_TEST_JULIA") else
    Sys.setenv(PDMPSAMPLERSR_TEST_JULIA = old), add = TRUE)
  Sys.setenv(PDMPSAMPLERSR_TEST_JULIA = "false")
  calls <- 0L
  local_mocked_bindings(julia_setup = function(...) {
    calls <<- calls + 1L
    stop("Julia setup must not be called")
  }, .package = "JuliaCall")
  probes <- new.env(parent = globalenv())
  sys.source(test_path("helper-availability.R"), envir = probes)
  expect_identical(calls, 0L)
  expect_false(probes$is_julia_available())
  expect_false(probes$is_pdmp_julia_backend_available())
  expect_identical(calls, 0L)
})

test_that("an unavailable Julia runtime is probed once without warnings", {
  old <- Sys.getenv("PDMPSAMPLERSR_TEST_JULIA", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("PDMPSAMPLERSR_TEST_JULIA") else
    Sys.setenv(PDMPSAMPLERSR_TEST_JULIA = old), add = TRUE)
  Sys.unsetenv("PDMPSAMPLERSR_TEST_JULIA")
  calls <- 0L
  local_mocked_bindings(julia_setup = function(...) {
    calls <<- calls + 1L
    warning("Julia dependency installation failed")
    stop("Julia is unavailable")
  }, .package = "JuliaCall")
  probes <- new.env(parent = globalenv())
  sys.source(test_path("helper-availability.R"), envir = probes)
  expect_warning(expect_false(probes$is_julia_available()), NA)
  expect_false(probes$is_pdmp_julia_backend_available())
  expect_false(probes$is_julia_available())
  expect_identical(calls, 1L)
})
