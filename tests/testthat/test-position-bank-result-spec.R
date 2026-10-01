test_that("supplied frozen bank survives unavailable live counter payload", {
  anchors <- matrix(seq_len(24), nrow = 3)
  supplied <- list(version = 1L, ordered_anchor_positions = anchors,
    distance_scales = c(1, 2, 3), selection_policy = "position",
    tie_policy = "lowest_index_outgoing", requested_capacity = 8L,
    parameter_order = letters[1:3], dataset_identity = list(hash = "data"),
    envelope_identity = list(hash = "envelope"),
    immutable_envelope_grouping_spec = list(old = TRUE))
  endpoint <- list(flow = list(mean = rep(0, 3), covariance = diag(3)))
  out <- .position_bank_result_spec(list(), supplied, endpoint, "position",
    letters[1:3], "data", "envelope", list(current = TRUE),
    function(...) stop("node support must not be requested"))
  expect_identical(out$ordered_anchor_positions, anchors)
  expect_identical(out$distance_scales, supplied$distance_scales)
  expect_identical(out$immutable_envelope_grouping_spec, list(current = TRUE))
  expect_identical(out$physical_endpoint, endpoint)
  expect_identical(out$frozen_flow, endpoint$flow)
})

test_that("returned live bank disagreement still fails closed", {
  supplied <- list(ordered_anchor_positions = matrix(1:6, nrow = 2),
    distance_scales = c(1, 2))
  expect_error(.position_bank_result_spec(
    list(anchor_positions = matrix(7:12, nrow = 2)), supplied, NULL,
    "position", letters[1:2], "data", "envelope", list(), identity),
    "live position bank differs")
  expect_error(.position_bank_result_spec(list(), NULL, NULL, "position",
    letters[1:2], "data", "envelope", list(), identity),
    "newly prepared position bank")
})

test_that("new HCV banks retain warmup-learning provenance", {
  counters <- list(
    anchor_positions = matrix(1:6, nrow = 2),
    anchor_distance_scales = c(1, 2),
    anchor_bank_requested_capacity = 3L,
    warmup_bank_learning_policy = "exact_variance",
    warmup_bank_candidate_count = 7L,
    warmup_bank_training_state_count = 7L,
    warmup_bank_objective_single = 12,
    warmup_bank_objective_selected = 3,
    warmup_bank_learning_seconds = 0.25)
  out <- .position_bank_result_spec(counters, NULL, NULL, "position",
    letters[1:2], "data", "envelope", list(), identity)
  expect_identical(out$warmup_bank_learning_policy, "exact_variance")
  expect_identical(out$warmup_bank_candidate_count, 7L)
  expect_identical(out$warmup_bank_training_state_count, 7L)
  expect_equal(out$warmup_bank_objective_single /
    out$warmup_bank_objective_selected, 4)
  expect_equal(out$warmup_bank_learning_seconds, 0.25)
  supplied <- out
  expect_error(.position_bank_result_spec(
    list(warmup_bank_learning_policy = "geometric"), supplied, NULL,
    "position", letters[1:2], "data", "envelope", list(), identity),
    "bank-learning policy differs")
})
