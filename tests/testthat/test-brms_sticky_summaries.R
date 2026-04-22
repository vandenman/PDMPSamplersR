# Unit tests for Phase 3 sticky summary helpers.
# Does not require Julia or a real PDMP fit.
# Uses a minimal mock brmsfit carrying a sticky attribute.

make_mock_sticky_fit <- function(
    incl_chain1,
    incl_chain2 = NULL,
    coef_names = paste0("x", seq_along(incl_chain1)),
    model_prior = PDMPSamplersR::bernoulli(0.5)
) {
    unc_names <- c("b.Intercept", paste0("b.", coef_names), "sigma")
    can_stick <- c(FALSE, rep(TRUE, length(incl_chain1)), FALSE)
    incl_list <- list(
        chain1 = setNames(incl_chain1, paste0("b.", coef_names))
    )
    if (!is.null(incl_chain2))
        incl_list$chain2 <- setNames(incl_chain2, paste0("b.", coef_names))
    fit <- structure(list(), class = c("sticky_brmsfit", "brmsfit"))
    # Expand scalar bernoulli prob to length d for consistency with brm_pdmp output
    d <- length(unc_names)
    if (inherits(model_prior, "bernoulli") && length(model_prior$prob) == 1L)
        model_prior <- PDMPSamplersR::bernoulli(prob = rep(model_prior$prob, d))
    attr(fit, "sticky") <- list(
        inclusion_probs = incl_list,
        can_stick       = can_stick,
        unc_names       = unc_names,
        model_prior     = model_prior
    )
    fit
}

# ──────────────────────────────────────────────────────────────────────────────
# inclusion_prob
# ──────────────────────────────────────────────────────────────────────────────

test_that("inclusion_prob.brmsfit returns named numeric in [0, 1] (single chain)", {
    fit <- make_mock_sticky_fit(c(0.9, 0.1, 0.8))
    ip  <- PDMPSamplersR::inclusion_prob(fit)
    expect_named(ip)
    expect_true(is.numeric(ip))
    expect_length(ip, 3L)
    expect_true(all(ip >= 0 & ip <= 1))
    expect_equal(names(ip), c("b.x1", "b.x2", "b.x3"))
})

test_that("inclusion_prob.brmsfit averages across chains", {
    fit <- make_mock_sticky_fit(c(0.8, 0.2), c(0.6, 0.4))
    ip  <- PDMPSamplersR::inclusion_prob(fit)
    expect_equal(ip, c(b.x1 = 0.7, b.x2 = 0.3))
})

test_that("inclusion_prob.brmsfit selects chain by index", {
    fit <- make_mock_sticky_fit(c(0.9, 0.1), c(0.5, 0.5))
    ip  <- PDMPSamplersR::inclusion_prob(fit, chain = 2L)
    expect_equal(ip, c(b.x1 = 0.5, b.x2 = 0.5))
})

test_that("inclusion_prob.brmsfit errors for missing chain", {
    fit <- make_mock_sticky_fit(c(0.9, 0.1))
    expect_error(PDMPSamplersR::inclusion_prob(fit, chain = 3L), "not found")
})

test_that("inclusion_prob errors on non-sticky brmsfit", {
    fit <- structure(list(), class = "brmsfit")
    expect_error(PDMPSamplersR::inclusion_prob(fit), "sticky")
})

test_that("inclusion_prob errors on non-brmsfit input", {
    expect_error(PDMPSamplersR::inclusion_prob(list()), class = "rlang_error")
})

# ──────────────────────────────────────────────────────────────────────────────
# median_probability_model
# ──────────────────────────────────────────────────────────────────────────────

test_that("median_probability_model returns names above threshold", {
    fit <- make_mock_sticky_fit(c(0.9, 0.3, 0.7))
    mpm <- PDMPSamplersR::median_probability_model(fit)
    expect_equal(mpm, c("b.x1", "b.x3"))
})

test_that("median_probability_model with custom threshold", {
    fit <- make_mock_sticky_fit(c(0.9, 0.3, 0.7))
    mpm <- PDMPSamplersR::median_probability_model(fit, threshold = 0.8)
    expect_equal(mpm, "b.x1")
})

test_that("median_probability_model with threshold = 0 returns all names", {
    fit <- make_mock_sticky_fit(c(0.9, 0.01))
    mpm <- PDMPSamplersR::median_probability_model(fit, threshold = 0)
    expect_equal(mpm, c("b.x1", "b.x2"))
})

test_that("median_probability_model errors on invalid threshold", {
    fit <- make_mock_sticky_fit(c(0.5, 0.5))
    expect_error(PDMPSamplersR::median_probability_model(fit, threshold = 1.5), "threshold")
    expect_error(PDMPSamplersR::median_probability_model(fit, threshold = "a"), "threshold")
})

# ──────────────────────────────────────────────────────────────────────────────
# plot_inclusion
# ──────────────────────────────────────────────────────────────────────────────

test_that("plot_inclusion returns a ggplot object", {
    skip_if_not_installed("ggplot2")
    fit <- make_mock_sticky_fit(c(0.9, 0.1, 0.7))
    expect_warning(
        p <- PDMPSamplersR::plot_inclusion(fit),
        NA  # no warning expected from plot_inclusion itself
    )
    expect_s3_class(p, "gg")
})

# ──────────────────────────────────────────────────────────────────────────────
# inclusion_bf
# ──────────────────────────────────────────────────────────────────────────────

test_that("inclusion_bf returns named numeric vector (bernoulli prior)", {
    fit <- make_mock_sticky_fit(c(0.9, 0.1, 0.5), model_prior = PDMPSamplersR::bernoulli(0.5))
    bf <- PDMPSamplersR::inclusion_bf(fit)
    expect_named(bf)
    expect_true(is.numeric(bf))
    expect_length(bf, 3L)
    expect_equal(names(bf), c("b.x1", "b.x2", "b.x3"))
    # With bernoulli(0.5): prior odds = 1, so BF = posterior odds
    expect_equal(bf[["b.x1"]], 0.9 / 0.1)
    expect_equal(bf[["b.x2"]], 0.1 / 0.9)
    expect_equal(bf[["b.x3"]], 1.0)
})

test_that("inclusion_bf returns correct values with non-uniform prior", {
    fit <- make_mock_sticky_fit(c(0.6), model_prior = PDMPSamplersR::bernoulli(0.3))
    bf <- PDMPSamplersR::inclusion_bf(fit)
    # BF = (0.6/0.4) / (0.3/0.7) = 1.5 / (3/7) = 1.5 * 7/3 = 3.5
    expect_equal(bf[["b.x1"]], (0.6 / 0.4) / (0.3 / 0.7))
})

test_that("inclusion_bf works with betabernoulli prior", {
    fit <- make_mock_sticky_fit(c(0.8, 0.2), model_prior = PDMPSamplersR::betabernoulli(1, 1))
    bf <- PDMPSamplersR::inclusion_bf(fit)
    # betabernoulli(1,1): pi_prior = 1/(1+1) = 0.5; prior odds = 1
    expect_equal(bf[["b.x1"]], 0.8 / 0.2)
    expect_equal(bf[["b.x2"]], 0.2 / 0.8)
})

test_that("inclusion_bf errors without model_prior in sticky metadata", {
    fit <- make_mock_sticky_fit(c(0.5, 0.5))
    attr(fit, "sticky")$model_prior <- NULL
    expect_error(PDMPSamplersR::inclusion_bf(fit), "model prior")
})

test_that("inclusion_bf errors on non-sticky fit", {
    fit <- structure(list(), class = "brmsfit")
    expect_error(PDMPSamplersR::inclusion_bf(fit), "sticky")
})

# ──────────────────────────────────────────────────────────────────────────────
# check_inclusion_stability
# ──────────────────────────────────────────────────────────────────────────────

test_that("check_inclusion_stability returns data frame with correct columns", {
    fit <- make_mock_sticky_fit(c(0.9, 0.1), c(0.7, 0.3))
    df  <- PDMPSamplersR::check_inclusion_stability(fit)
    expect_s3_class(df, "data.frame")
    expect_named(df, c("coefficient", "mean", "sd", "min", "max", "range"))
    expect_equal(nrow(df), 2L)
})

test_that("check_inclusion_stability gives correct cross-chain statistics", {
    fit <- make_mock_sticky_fit(c(0.8, 0.2), c(0.6, 0.4))
    df  <- PDMPSamplersR::check_inclusion_stability(fit)
    expect_equal(df$mean, c(0.7, 0.3))
    expect_equal(df$range, c(0.2, 0.2))
})

test_that("check_inclusion_stability informs on single chain", {
    fit <- make_mock_sticky_fit(c(0.9, 0.1))
    expect_message(PDMPSamplersR::check_inclusion_stability(fit), "one chain")
})

test_that("check_inclusion_stability errors on non-sticky fit", {
    fit <- structure(list(), class = "brmsfit")
    expect_error(PDMPSamplersR::check_inclusion_stability(fit), "sticky")
})
