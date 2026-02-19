test_that("write_stan_json writes valid JSON", {
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file))

  data <- list(N = 10L, x = c(1.0, 2.0, 3.0))
  write_stan_json(data, file)

  expect_true(file.exists(file))
  parsed <- jsonlite::fromJSON(file)
  expect_equal(parsed$N, 10)
  expect_equal(parsed$x, c(1, 2, 3))
})

test_that("write_stan_json handles matrices", {
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file))

  m <- matrix(1:6, nrow = 2, ncol = 3)
  write_stan_json(list(m = m), file)

  parsed <- jsonlite::fromJSON(file)
  expect_equal(dim(parsed$m), c(2, 3))
})

test_that("write_stan_json converts logicals to integers", {
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file))

  write_stan_json(list(z = c(TRUE, FALSE, TRUE)), file)
  parsed <- jsonlite::fromJSON(file)
  expect_equal(parsed$z, c(1, 0, 1))
})

test_that("write_stan_json rejects invalid inputs", {
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file))

  expect_error(write_stan_json("not a list", file), "list")
  expect_error(write_stan_json(list(1, 2), file), "names")
  expect_error(write_stan_json(list(a = 1, a = 2), file), "Duplicate")
  expect_error(write_stan_json(list(a = NA_real_), file), "NA")
})
