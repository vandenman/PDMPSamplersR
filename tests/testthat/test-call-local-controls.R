test_that("call-local control constructors validate settings", {
  expect_identical(sampler_control()$hcv_cache, TRUE)
  expect_null(sampler_control()$warmup_stop)
  expect_identical(sampler_control(warmup_stop = list())$warmup_stop, list())
  expect_identical(sampler_control(grid_n_min = 3)$grid_n_min, 3L)
  expect_error(sampler_control(grid_n_min = 0), "positive integer")
  expect_error(sampler_control(profile = NA), "TRUE or FALSE")
  expect_error(sampler_control(warmup_stop = list(check_every = 0)), "Invalid warmup")
  expect_error(sampler_control(warmup_stop = list(typo = 1)), "named list")
  expect_error(sampler_control(warmup_stop = list(1)), "named list")
  expect_error(stan_compile_options(native = NA), "TRUE or FALSE")
  expect_error(stan_compile_options(cache_dir = ""), "empty")
  expect_identical(stan_compile_options(cache = TRUE)$cache, TRUE)
})

test_that("project selection uses the R option instead of ambient environment", {
  old <- Sys.getenv("PDMPSAMPLERSR_CONFIG_JULIA_PROJECT_DIR", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("PDMPSAMPLERSR_CONFIG_JULIA_PROJECT_DIR") else
    Sys.setenv(PDMPSAMPLERSR_CONFIG_JULIA_PROJECT_DIR = old))
  Sys.setenv(PDMPSAMPLERSR_CONFIG_JULIA_PROJECT_DIR = "must-not-select-this-project")
  previous <- options(pdmpsamplersr_config_julia_project_dir = NULL)
  on.exit(options(previous), add = TRUE)
  expect_identical(get_julia_project_dir(), default_julia_project_dir())
  options(pdmpsamplersr_config_julia_project_dir = "/explicit/project")
  expect_identical(get_julia_project_dir(), "/explicit/project")
})

test_that("online summary options do not leak between calls in one worker", {
  skip_if_not(is_pdmp_julia_backend_available(), "Julia backend unavailable")
  variables <- c("PDMPSAMPLERS_ONLINE_BATCHES", "PDMPSAMPLERS_ONLINE_GRID_SPACING",
    "PDMPSAMPLERS_STREAM_EVENTS", "PDMPSAMPLERS_ONLINE_WARMUP_BATCHES",
    "PDMP_GRID_N_MIN", "PDMP_ADAPTIVE_BOOMERANG_LAMBDA_REF")
  old <- Sys.getenv(variables, unset = NA_character_)
  on.exit(for (name in variables) {
    if (is.na(old[[name]])) Sys.unsetenv(name) else
      do.call(Sys.setenv, setNames(list(old[[name]]), name))
  }, add = TRUE)
  do.call(Sys.setenv, setNames(as.list(rep("invalid-ambient-setting", length(variables))), variables))
  run <- function(path, online) pdmp_sample(function(x) x, d = 2,
    flow = "BouncyParticle", algorithm = "GridThinningStrategy", T = 0.4,
    t_warmup = 0.1, seed = 123, x0 = c(0.2, -0.1),
    trace_storage = "streaming", streaming_trace_directory = path,
    online_summaries = online, show_progress = FALSE)
  first_dir <- tempfile(); second_dir <- tempfile(); third_dir <- tempfile()
  on.exit(unlink(c(first_dir, second_dir, third_dir), recursive = TRUE), add = TRUE)
  first <- run(first_dir, online_summary_options(batches = 4L,
    grid_spacing = 0.05, warmup_batches = 2L))
  first_summary <- read_online_summaries(first_dir)
  expect_equal(first_summary$n_batches, 4)
  expect_equal(read_online_summaries(first_dir, "warmup")$n_batches, 2)
  expect_gt(length(first_summary$grid_times), 0)
  second <- run(second_dir, NULL)
  expect_length(list.files(second_dir, "online_summaries.bin", recursive = TRUE), 0L)
  third <- run(third_dir, online_summary_options(batches = 2L, write_events = TRUE))
  expect_equal(read_online_summaries(third_dir)$n_batches, 2)
  expect_identical(read_online_summaries(first_dir), first_summary)
  expect_true(all(Sys.getenv(variables) == "invalid-ambient-setting"))
  cleaned <- cleanup_worker_state()
  expect_true(unname(as.logical(cleaned)))
  expect_identical(attr(cleaned, "reason"), "cleaned")
  expect_identical(read_online_summaries(first_dir), first_summary)
})

test_that("Stan profiling and streaming controls apply to one call", {
  skip_if_not(is_pdmp_julia_backend_available(), "Julia backend unavailable")
  model <- system.file("stan", "models", "mvnormal.stan", package = "PDMPSamplersR")
  data <- list(N = 2L, mu = c(0, 0), sigma = diag(2))
  directory <- tempfile()
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  args <- list(path_to_stanmodel = model, standata = data,
    flow = "BouncyParticle", algorithm = "GridThinningStrategy", T = 0.4,
    t_warmup = 0.1, x0 = c(0.2, -0.1), seed = 44L, show_progress = FALSE)
  fit <- do.call(pdmp_sample_from_stanmodel, c(args, list(
    control = sampler_control(profile = TRUE, grid_n_min = 2L),
    trace_storage = "streaming", streaming_trace_directory = directory,
    online_summaries = online_summary_options(batches = 2L))))
  expect_equal(read_online_summaries(directory)$n_batches, 2)
  profile <- attr(fit, "bridge_phase_profile")
  expect_true(is.list(profile))
  expect_true(all(c("julia_construction_seconds", "julia_sample_seconds") %in%
                    names(profile$julia_internal)))
  ordinary <- do.call(pdmp_sample_from_stanmodel, args)
  expect_null(attr(ordinary, "bridge_phase_profile"))
})
