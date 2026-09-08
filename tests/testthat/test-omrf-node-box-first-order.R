test_that("node-box first-order specifications stay compact and deterministic", {
  set.seed(41)
  X <- matrix(sample.int(4L, 500L * 30L, replace = TRUE) - 1L,
              nrow = 500L, ncol = 30L)
  build <- function(groups, prior, geometry = "node_box_first_order") omrf_residual_envelope(
    X, rep.int(4L, 30L), "thresholds_0", "interactions_0",
    backend = "analytic", use_hcv = FALSE,
    factorization = "person_node", factor_sampling = "size_biased",
    prior_backend = prior, geometry = geometry,
    box_groups = groups)
  for (prior in c("stan", "analytic")) {
    envelope <- build(4L, prior)
    expect_identical(envelope$bound_type, "node_box_first_order")
    expect_null(envelope$weights)
    expect_null(envelope$node_weights)
    expect_identical(dim(envelope$node_box_membership), c(30L, 500L))
    expect_identical(ncol(envelope$node_box_minima), 30L)
    expect_identical(dim(envelope$node_box_minima),
                     dim(envelope$node_box_maxima))
    expect_true(isTRUE(envelope$exact_selected_refinement))
  }
  first <- build(8L, "analytic")
  second <- build(8L, "analytic")
  expect_identical(first$node_box_membership, second$node_box_membership)
  expect_identical(first$node_box_minima, second$node_box_minima)
  expect_identical(first$node_box_maxima, second$node_box_maxima)
  # Membership plus extrema is well below the forbidden dense P-by-NP matrix.
  expect_lt(as.numeric(object.size(first)), 50 * 1024^2)

  cell_marked <- build(2L, "analytic",
                       "cell_marked_node_box_first_order")
  expect_identical(cell_marked$bound_type,
                   "cell_marked_node_box_first_order")
  expect_null(cell_marked$weights)
  expect_null(cell_marked$node_weights)
  expect_true(isTRUE(cell_marked$exact_selected_refinement))
  expect_identical(length(cell_marked$node_box_nodes), 60L)
  expect_error(build(4L, "analytic",
                     "cell_marked_node_box_first_order"),
               "requires.*box_groups = 2")

  malformed <- build(8L, "stan")
  malformed$node_box_membership <- malformed$node_box_membership[-1, ]
  expect_error(stan_subsampling(500L * 30L, 1L,
    prior_standata = list(dummy = 1), residual_envelope = malformed),
    "P-by-N integer membership")
  expect_error(stan_subsampling(500L * 30L + 1L, 1L,
    prior_standata = list(dummy = 1), residual_envelope = build(2L, "stan")),
    "requires 15000 likelihood contributions")
})
