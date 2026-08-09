benchmark_script <- function(name) {
  source_path <- testthat::test_path("..", "..", "benchmarks", name)
  if (file.exists(source_path)) return(source_path)
  installed_path <- system.file("benchmarks", name, package = "PDMPSamplersR")
  if (!nzchar(installed_path)) stop("Installed benchmark script not found: ", name)
  installed_path
}

benchmark_implementation <- function(name) {
  source_path <- testthat::test_path("..", "..", "inst", "benchmarks", name)
  if (file.exists(source_path)) return(source_path)
  installed_path <- system.file("benchmarks", name, package = "PDMPSamplersR")
  if (!nzchar(installed_path)) stop("Benchmark implementation not found: ", name)
  installed_path
}

test_that("subsampling benchmark summaries are paired and report uncertainty", {
  environment <- new.env(parent = baseenv())
  sys.source(
    benchmark_script("summarize_subsampling_first_batch.R"),
    envir = environment
  )
  rows <- expand.grid(
    method = c("full", "subsampling"), seed = 1:3,
    stringsAsFactors = FALSE
  )
  rows$family <- "binomial"
  rows$flow <- "ZigZag"
  rows$sticky <- FALSE
  rows$anchor_source <- "estimated"
  rows$N <- 500L
  rows$dimension <- 3L
  rows$comparison_minibatch <- 50L
  rows$benchmark_mode <- "adaptive_diagnostics"
  rows$pdmp_duration <- ifelse(rows$method == "subsampling", 200, 100)
  rows$sampling_seconds <- ifelse(rows$method == "subsampling", 2, 4)
  rows$wall_seconds <- ifelse(rows$method == "subsampling", 5, 8)
  rows$min_ess_per_sampling_second <- ifelse(rows$method == "subsampling", 20, 10)
  rows$min_ess_per_wall_second <- ifelse(rows$method == "subsampling", 8, 5)
  rows$ess_diagnostic_reliable <- TRUE
  rows$candidates_per_event <- ifelse(rows$method == "subsampling", 4, NA)
  rows$cell_roof_proposals_per_event <- ifelse(rows$method == "subsampling", 5, NA)
  rows$aggregate_acceptance_rate <- ifelse(rows$method == "subsampling", 0.8, NA)
  rows$final_acceptance_rate <- ifelse(rows$method == "subsampling", 0.25, NA)
  rows$subset_evaluations <- ifelse(rows$method == "subsampling", 40, NA)

  result <- environment$summarize_subsampling_first_batch(rows)
  expect_equal(nrow(result$summary), 1L)
  expect_equal(result$summary$n_pairs, 3L)
  expect_equal(result$summary$n_ess_pairs, 3L)
  expect_equal(result$summary$sampling_pdmp_throughput_full_median, 25)
  expect_equal(result$summary$sampling_pdmp_throughput_subsampling_median, 100)
  expect_equal(result$summary$sampling_pdmp_throughput_ratio_median, 4)
  expect_equal(result$summary$wall_pdmp_throughput_ratio_median, 3.2)
  expect_true(all(c("sampling_pdmp_throughput_ratio_q25",
                    "sampling_pdmp_throughput_ratio_q75") %in%
                    names(result$summary)))
})

test_that("multiple minibatches pair one-to-one and remain separate", {
  environment <- new.env(parent = baseenv())
  sys.source(
    benchmark_script("summarize_subsampling_first_batch.R"),
    envir = environment
  )
  rows <- expand.grid(
    method = c("full", "subsampling"), comparison_minibatch = c(25L, 50L),
    stringsAsFactors = FALSE
  )
  rows$family <- "binomial"
  rows$flow <- "ZigZag"
  rows$sticky <- FALSE
  rows$anchor_source <- "estimated"
  rows$seed <- 1L
  rows$N <- 500L
  rows$dimension <- 3L
  rows$benchmark_mode <- "equal_duration"
  rows$pdmp_duration <- 100
  rows$sampling_seconds <- ifelse(rows$method == "subsampling", 2, 4)
  rows$wall_seconds <- ifelse(rows$method == "subsampling", 5, 8)
  rows$min_ess_per_sampling_second <- ifelse(rows$method == "subsampling", 20, 10)
  rows$min_ess_per_wall_second <- ifelse(rows$method == "subsampling", 8, 5)
  rows$ess_diagnostic_reliable <- TRUE
  rows$candidates_per_event <- ifelse(rows$method == "subsampling", 4, NA)
  rows$cell_roof_proposals_per_event <- ifelse(rows$method == "subsampling", 5, NA)
  rows$aggregate_acceptance_rate <- ifelse(rows$method == "subsampling", 0.8, NA)
  rows$final_acceptance_rate <- ifelse(rows$method == "subsampling", 0.25, NA)
  rows$subset_evaluations <- ifelse(rows$method == "subsampling", 40, NA)

  result <- environment$summarize_subsampling_first_batch(rows)
  expect_equal(nrow(result$paired), 2L)
  expect_equal(nrow(result$summary), 2L)
  expect_equal(sort(result$summary$comparison_minibatch), c(25L, 50L))
  expect_equal(result$summary$n_pairs, c(1L, 1L))
  expect_equal(result$summary$equal_duration_wall_speedup_median, c(1.6, 1.6))

  duplicated <- rbind(rows, rows[rows$method == "full" &
                                  rows$comparison_minibatch == 25L, ])
  expect_error(
    environment$summarize_subsampling_first_batch(duplicated),
    "exactly one full row"
  )

  missing_subsampling <- rows[!(rows$method == "subsampling" &
                            rows$comparison_minibatch == 50L), ]
  expect_error(
    environment$summarize_subsampling_first_batch(missing_subsampling),
    "one-to-one pairing"
  )
})

test_that("HCV benchmark summaries use reliable paired triplets", {
  environment <- new.env(parent = baseenv())
  sys.source(
    benchmark_script("summarize_subsampling_first_batch.R"),
    envir = environment
  )
  sys.source(
    benchmark_script("summarize_subsampling_hcv.R"),
    envir = environment
  )
  rows <- expand.grid(
    method = c("full", "subsampling", "hcv"), seed = 1:3,
    stringsAsFactors = FALSE
  )
  rows$family <- "bernoulli"; rows$flow <- "BouncyParticle"
  rows$sticky <- FALSE; rows$anchor_source <- "estimated"
  rows$N <- 500L; rows$dimension <- 6L; rows$comparison_minibatch <- 50L
  rows$benchmark_mode <- "equal_duration"; rows$pdmp_duration <- 80
  rows$sampling_seconds <- 1; rows$wall_seconds <- 2
  subsampling_ess <- c(1, 50, 100); hcv_ess <- c(10, 20, 1000)
  rows$min_ess_per_sampling_second <- 1
  rows$min_ess_per_wall_second <- 0.5
  for (seed in 1:3) {
    rows$min_ess_per_sampling_second[rows$method == "subsampling" & rows$seed == seed] <-
      subsampling_ess[[seed]]
    rows$min_ess_per_sampling_second[rows$method == "hcv" & rows$seed == seed] <-
      hcv_ess[[seed]]
  }
  rows$min_ess_per_wall_second <- rows$min_ess_per_sampling_second / 2
  rows$ess_diagnostic_reliable <- TRUE
  rows$ess_diagnostic_reliable[rows$method == "hcv" & rows$seed == 2] <- FALSE
  rows$cell_roof_proposals_per_event <- ifelse(rows$method == "full", NA, 4)
  rows$final_acceptance_rate <- ifelse(rows$method == "full", NA, 0.25)
  rows$subset_evaluations <- ifelse(rows$method == "full", NA, 40)
  rows$grid_bound_violations <- 0

  result <- environment$summarize_subsampling_hcv(rows)
  expect_equal(nrow(result$paired), 3L)
  expect_equal(result$summary$n_triplets, 3L)
  expect_equal(result$summary$n_ess_pairs_hcv_vs_subsampling, 2L)
  expect_true(is.na(result$paired$min_ess_sampling_ratio_hcv_vs_subsampling[[2L]]))
  expect_equal(result$summary$min_ess_sampling_ratio_hcv_vs_subsampling_median, 10)
  expect_false(isTRUE(all.equal(
    result$summary$min_ess_sampling_ratio_hcv_vs_subsampling_median,
    stats::median(hcv_ess) / stats::median(subsampling_ess)
  )))

  duplicated <- rbind(rows, rows[rows$method == "hcv" & rows$seed == 1, ])
  expect_error(environment$summarize_subsampling_hcv(duplicated),
               "exactly one hcv row")
  missing <- rows[!(rows$method == "subsampling" & rows$seed == 3), ]
  expect_error(environment$summarize_subsampling_hcv(missing),
               "one-to-one pairing")
})

test_that("OMRF rows pair current subsampling method and summarize proposal diagnostics", {
  environment <- new.env(parent = baseenv())
  sys.source(
    benchmark_script("summarize_subsampling_first_batch.R"),
    envir = environment)
  rows <- expand.grid(
    method = c("full", "subsampling"), N = c(200L, 500L),
    P = c(10L, 20L, 30L), seed = 1L,
    stringsAsFactors = FALSE)
  rows$dimension <- 2L * rows$P + rows$P * (rows$P - 1L) / 2L
  rows$comparison_minibatch <- ifelse(rows$N == 200L, 20L, 50L)
  rows$benchmark_mode <- "calibrated_equal_duration"
  rows$pdmp_duration <- 0.01
  rows$sampling_seconds <- ifelse(rows$method == "subsampling", 2, 1)
  rows$wall_seconds <- ifelse(rows$method == "subsampling", 3, 2)
  rows$min_ess_per_sampling_second <- NA_real_
  rows$min_ess_per_wall_second <- NA_real_
  rows$ess_diagnostic_reliable <- FALSE
  rows$subsampling_cell_roof_proposals <- ifelse(rows$method == "subsampling", 100, NA)
  rows$subsampling_aggregate_accepts <- ifelse(rows$method == "subsampling", 40, NA)
  rows$subsampling_subset_evaluations <- ifelse(rows$method == "subsampling", 40, NA)
  rows$subsampling_final_reflections <- ifelse(rows$method == "subsampling", 10, NA)
  rows$selected_gradient_calls <- ifelse(rows$method == "subsampling", 80, NA)
  rows$persons_evaluated <- ifelse(
    rows$method == "subsampling",
    80 * rows$comparison_minibatch, NA)

  result <- environment$summarize_subsampling_first_batch(rows)
  expect_equal(nrow(result$paired), 6L)
  expect_equal(nrow(result$summary), 6L)
  expect_equal(sort(unique(result$summary$P)), c(10L, 20L, 30L))
  expect_equal(result$paired$roof_proposals_per_subset, rep(2.5, 6L))
  expect_equal(result$paired$aggregate_acceptance_rate, rep(0.4, 6L))
  expect_equal(result$paired$final_acceptance_rate, rep(0.25, 6L))
  expect_equal(
    result$paired$persons_per_selected_gradient,
    result$paired$comparison_minibatch)
  expect_equal(result$summary$n_ess_pairs, rep(0L, 6L))
})

test_that("OMRF preflight projects reliable horizons and complete paired cost", {
  environment <- new.env(parent = globalenv())
  sys.source(
    benchmark_implementation("omrf_subsampling_impl.R"), envir = environment)
  probes <- expand.grid(
    method = c("full", "subsampling"), seed = 1:3,
    pdmp_duration = c(0.25, 0.5), stringsAsFactors = FALSE)
  probes$N <- 200L
  probes$P <- 10L
  probes$comparison_minibatch <- 20L
  probes$selected_gradient_calls <- ifelse(
    probes$method == "subsampling",
    ifelse(probes$pdmp_duration == 0.25, 1000, 4000), NA_real_)
  probes$persons_evaluated <-
    probes$selected_gradient_calls * probes$comparison_minibatch
  probes$main_events <- ifelse(
    probes$method == "subsampling",
    ifelse(probes$pdmp_duration == 0.25, 100, 200),
    ifelse(probes$pdmp_duration == 0.25, 150, 300))
  probes$initialization_seconds <- ifelse(probes$method == "subsampling", 0.2, 0.3)
  probes$sampling_seconds <- ifelse(
    probes$method == "subsampling",
    ifelse(probes$pdmp_duration == 0.25, 1, 4),
    ifelse(probes$pdmp_duration == 0.25, 0.5, 1))
  probes$wall_seconds <- probes$initialization_seconds + probes$sampling_seconds
  probes$anchor_preparation_seconds <- 0.05

  projection <- environment$project_omrf_preflight(
    probes, maximum_grid_wall_seconds = 1e6)
  expect_equal(nrow(projection), 1L)
  expect_equal(projection$replicate_count, 3L)
  expect_equal(projection$projected_selected_gradient_calls, 200000,
               tolerance = 1e-8)
  expect_equal(
    projection$projected_person_evaluations,
    20 * projection$projected_selected_gradient_calls)
  expect_true(projection$calibrated_physical_duration >= 1)
  expect_true(projection$projected_subsampling_main_events >= 500)
  expect_true(projection$projected_full_main_events >= 500)
  expect_true(projection$projection_supported)
  expect_true(projection$ess_reliability_eligible)
  expect_true(projection$preflight_feasible)
  expect_equal(
    projection$projected_complete_grid_wall_seconds,
    3 * projection$projected_pair_wall_seconds)
  expect_false(projection$ess_projection_available)

  infeasible <- environment$project_omrf_preflight(
    probes, maximum_duration = 2, maximum_grid_wall_seconds = 1e6)
  expect_false(infeasible$within_maximum_duration)
  expect_false(infeasible$preflight_feasible)
  expect_error(
    environment$assert_omrf_preflight_feasible(
      infeasible, tempfile(fileext = ".csv")),
    "No acceptance sampling was started")
})

test_that("benchmark-only OMRF anchor gradient matches finite differences", {
  environment <- new.env(parent = globalenv())
  sys.source(
    benchmark_implementation("omrf_subsampling_impl.R"), envir = environment)
  stan_data <- environment$simulate_omrf_benchmark_data(20L, 4L, 91L)
  q <- seq(-0.2, 0.2, length.out = environment$.omrf_dimension(4L))
  analytic <- environment$.omrf_log_posterior(q, stan_data, gradient = TRUE)
  epsilon <- 1e-6
  numerical <- vapply(seq_along(q), function(index) {
    plus <- minus <- q
    plus[index] <- plus[index] + epsilon
    minus[index] <- minus[index] - epsilon
    (environment$.omrf_log_posterior(plus, stan_data) -
       environment$.omrf_log_posterior(minus, stan_data)) / (2 * epsilon)
  }, numeric(1L))
  expect_equal(analytic, numerical, tolerance = 2e-6)
})

test_that("benchmark wrappers are synchronized and implementations are canonical", {
  source_root <- testthat::test_path("..", "..")
  if (!file.exists(file.path(source_root, "benchmarks", "omrf_subsampling.R"))) {
    source_root <- file.path(source_root, "00_pkg_src", "PDMPSamplersR")
  }
  installed_root <- system.file("benchmarks", package = "PDMPSamplersR")
  for (name in c("omrf_subsampling.R", "summarize_subsampling_first_batch.R",
                 "summarize_subsampling_hcv.R")) {
    expect_identical(
      readLines(file.path(source_root, "benchmarks", name), warn = FALSE),
      readLines(file.path(installed_root, name), warn = FALSE))
  }
  expect_true(file.exists(file.path(installed_root, "omrf_subsampling_impl.R")))
  expect_true(file.exists(file.path(
    installed_root, "summarize_subsampling_first_batch_impl.R")))
  expect_true(file.exists(file.path(
    installed_root, "summarize_subsampling_hcv_impl.R")))
})
