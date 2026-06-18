
rethrow_pdmp_julia_error <- function(err) {
  msg <- conditionMessage(err)
  if (grepl("SupportBoundaryError:", msg, fixed = TRUE)) {
    stop(structure(
      list(message = msg, call = conditionCall(err), parent = err),
      class = c("pdmp_support_boundary_error", class(err))
    ))
  }
  stop(err)
}

.pdmpsamplers_julia_eval <- function(...) {
  # a very simple wrapper around JuliaCall for two reasons:

  # 1.
  # on MacOS CI, testthat testthat::local_mocked_bindings
  # segfaults when directly mocking JuliaCall::julia_call
  # this allows us to mock this wrapper instead.

  # 2.
  # this allows us to catch error generically and rethrow with a custom class for support boundary errors

  tryCatch(
    JuliaCall::julia_eval(...),
    error = rethrow_pdmp_julia_error
  )
}

.pdmpsamplers_julia_call <- function(...) {
  # a very simple wrapper around JuliaCall for two reasons:

  # 1.
  # on MacOS CI, testthat testthat::local_mocked_bindings
  # segfaults when directly mocking JuliaCall::julia_call
  # this allows us to mock this wrapper instead.
  JuliaCall::julia_call(...)
}
