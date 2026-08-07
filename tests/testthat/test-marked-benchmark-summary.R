test_that("marked benchmark summaries are paired and report uncertainty", {
  environment <- new.env(parent = baseenv())
  sys.source(
    testthat::test_path("..", "..", "benchmarks", "summarize_marked_first_batch.R"),
    envir = environment
  )
  rows <- expand.grid(
    method = c("full", "marked"), seed = 1:3,
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
  rows$pdmp_duration <- ifelse(rows$method == "marked", 200, 100)
  rows$sampling_seconds <- ifelse(rows$method == "marked", 2, 4)
  rows$wall_seconds <- ifelse(rows$method == "marked", 5, 8)
  rows$min_ess_per_sampling_second <- ifelse(rows$method == "marked", 20, 10)
  rows$min_ess_per_wall_second <- ifelse(rows$method == "marked", 8, 5)
  rows$ess_diagnostic_reliable <- TRUE
  rows$candidates_per_event <- ifelse(rows$method == "marked", 4, NA)
  rows$cell_roof_proposals_per_event <- ifelse(rows$method == "marked", 5, NA)
  rows$aggregate_acceptance_rate <- ifelse(rows$method == "marked", 0.8, NA)
  rows$final_acceptance_rate <- ifelse(rows$method == "marked", 0.25, NA)
  rows$subset_evaluations <- ifelse(rows$method == "marked", 40, NA)

  result <- environment$summarize_marked_first_batch(rows)
  expect_equal(nrow(result$summary), 1L)
  expect_equal(result$summary$n_pairs, 3L)
  expect_equal(result$summary$n_ess_pairs, 3L)
  expect_equal(result$summary$sampling_pdmp_throughput_full_median, 25)
  expect_equal(result$summary$sampling_pdmp_throughput_marked_median, 100)
  expect_equal(result$summary$sampling_pdmp_throughput_ratio_median, 4)
  expect_equal(result$summary$wall_pdmp_throughput_ratio_median, 3.2)
  expect_true(all(c("sampling_pdmp_throughput_ratio_q25",
                    "sampling_pdmp_throughput_ratio_q75") %in%
                    names(result$summary)))
})

test_that("multiple minibatches pair one-to-one and remain separate", {
  environment <- new.env(parent = baseenv())
  sys.source(
    testthat::test_path("..", "..", "benchmarks", "summarize_marked_first_batch.R"),
    envir = environment
  )
  rows <- expand.grid(
    method = c("full", "marked"), comparison_minibatch = c(25L, 50L),
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
  rows$sampling_seconds <- ifelse(rows$method == "marked", 2, 4)
  rows$wall_seconds <- ifelse(rows$method == "marked", 5, 8)
  rows$min_ess_per_sampling_second <- ifelse(rows$method == "marked", 20, 10)
  rows$min_ess_per_wall_second <- ifelse(rows$method == "marked", 8, 5)
  rows$ess_diagnostic_reliable <- TRUE
  rows$candidates_per_event <- ifelse(rows$method == "marked", 4, NA)
  rows$cell_roof_proposals_per_event <- ifelse(rows$method == "marked", 5, NA)
  rows$aggregate_acceptance_rate <- ifelse(rows$method == "marked", 0.8, NA)
  rows$final_acceptance_rate <- ifelse(rows$method == "marked", 0.25, NA)
  rows$subset_evaluations <- ifelse(rows$method == "marked", 40, NA)

  result <- environment$summarize_marked_first_batch(rows)
  expect_equal(nrow(result$paired), 2L)
  expect_equal(nrow(result$summary), 2L)
  expect_equal(sort(result$summary$comparison_minibatch), c(25L, 50L))
  expect_equal(result$summary$n_pairs, c(1L, 1L))
  expect_equal(result$summary$equal_duration_wall_speedup_median, c(1.6, 1.6))

  duplicated <- rbind(rows, rows[rows$method == "full" &
                                  rows$comparison_minibatch == 25L, ])
  expect_error(
    environment$summarize_marked_first_batch(duplicated),
    "exactly one full row"
  )

  missing_marked <- rows[!(rows$method == "marked" &
                            rows$comparison_minibatch == 50L), ]
  expect_error(
    environment$summarize_marked_first_batch(missing_marked),
    "one-to-one pairing"
  )
})

test_that("HCV benchmark summaries use reliable paired triplets", {
  environment <- new.env(parent = baseenv())
  sys.source(
    testthat::test_path("..", "..", "benchmarks", "summarize_marked_first_batch.R"),
    envir = environment
  )
  sys.source(
    testthat::test_path("..", "..", "benchmarks", "summarize_marked_hcv.R"),
    envir = environment
  )
  rows <- expand.grid(
    method = c("full", "marked", "hcv"), seed = 1:3,
    stringsAsFactors = FALSE
  )
  rows$family <- "bernoulli"; rows$flow <- "BouncyParticle"
  rows$sticky <- FALSE; rows$anchor_source <- "estimated"
  rows$N <- 500L; rows$dimension <- 6L; rows$comparison_minibatch <- 50L
  rows$benchmark_mode <- "equal_duration"; rows$pdmp_duration <- 80
  rows$sampling_seconds <- 1; rows$wall_seconds <- 2
  marked_ess <- c(1, 50, 100); hcv_ess <- c(10, 20, 1000)
  rows$min_ess_per_sampling_second <- 1
  rows$min_ess_per_wall_second <- 0.5
  for (seed in 1:3) {
    rows$min_ess_per_sampling_second[rows$method == "marked" & rows$seed == seed] <-
      marked_ess[[seed]]
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

  result <- environment$summarize_marked_hcv(rows)
  expect_equal(nrow(result$paired), 3L)
  expect_equal(result$summary$n_triplets, 3L)
  expect_equal(result$summary$n_ess_pairs_hcv_vs_marked, 2L)
  expect_true(is.na(result$paired$min_ess_sampling_ratio_hcv_vs_marked[[2L]]))
  expect_equal(result$summary$min_ess_sampling_ratio_hcv_vs_marked_median, 10)
  expect_false(isTRUE(all.equal(
    result$summary$min_ess_sampling_ratio_hcv_vs_marked_median,
    stats::median(hcv_ess) / stats::median(marked_ess)
  )))

  duplicated <- rbind(rows, rows[rows$method == "hcv" & rows$seed == 1, ])
  expect_error(environment$summarize_marked_hcv(duplicated),
               "exactly one hcv row")
  missing <- rows[!(rows$method == "marked" & rows$seed == 3), ]
  expect_error(environment$summarize_marked_hcv(missing),
               "one-to-one pairing")
})
