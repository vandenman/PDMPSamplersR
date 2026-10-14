test_that("supplied frozen bank survives unavailable live counter payload", {
  anchors <- matrix(seq_len(24), nrow = 3)
  supplied <- list(version = 1L, ordered_anchor_positions = anchors,
    distance_scales = c(1, 2, 3), selection_policy = "position",
    tie_policy = "lowest_index_outgoing", requested_capacity = 8L,
    parameter_order = letters[1:3], dataset_identity = list(hash = "data"),
    envelope_identity = list(hash = "envelope"),
    immutable_envelope_grouping_spec = list(old = TRUE))
  endpoint <- list(flow = list(mean = rep(0, 3), covariance = diag(3)))
  out <- .position_bank_result_spec(list(), supplied, endpoint,
    letters[1:3], "data", "envelope", list(current = TRUE))
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
    letters[1:2], "data", "envelope", list()),
    "live position bank differs")
  expect_error(.position_bank_result_spec(list(), NULL, NULL,
    letters[1:2], "data", "envelope", list()),
    "newly prepared position bank")
})
