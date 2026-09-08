test_that("compact OMRF norm specifications scale with NP", {
  make_x <- function(N, P) matrix(rep(0:2, length.out=N*P), N, P)
  build <- function(N, P, hcv=FALSE) omrf_residual_envelope(
    make_x(N,P), rep.int(3L,P), "thresholds_0", "interactions_0",
    backend="analytic", use_hcv=hcv,
    hcv_damping=if(hcv) Inf else 10, factorization="person_node",
    factor_sampling="size_biased", prior_backend="analytic",
    geometry=if(hcv) "norm_local_hcv_compact" else "norm_local")
  small<-build(100L,10L); large<-build(200L,10L); hcv<-build(100L,10L,TRUE)
  expect_null(small$weights); expect_null(small$node_weights)
  expect_null(hcv$weights); expect_null(hcv$node_weights)
  expect_false(any(c("weights", "node_weights") %in% names(small)))
  expect_false(any(c("weights", "node_weights") %in% names(hcv)))
  expect_equal(dim(small$norm_weights),c(30L,100L))
  expect_equal(dim(hcv$norm_weights),c(40L,100L))
  expect_lt(as.numeric(object.size(large))/as.numeric(object.size(small)),2.2)
  expect_gt(as.numeric(object.size(large))/as.numeric(object.size(small)),1.7)
})
