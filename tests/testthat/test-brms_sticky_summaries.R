# Unit tests for Phase 3 sticky summary helpers.
# Does not require Julia or a real PDMP fit.
# Uses a minimal mock brmsfit carrying a sticky attribute.

make_mock_sticky_fit <- function(
    incl_chain1,
    incl_chain2 = NULL,
    coef_names = paste0("x", seq_along(incl_chain1))
) {
    unc_names <- c("b.Intercept", paste0("b.", coef_names), "sigma")
    can_stick <- c(FALSE, rep(TRUE, length(incl_chain1)), FALSE)
    incl_list <- list(
        chain1 = setNames(incl_chain1, paste0("b.", coef_names))
    )
    if (!is.null(incl_chain2))
        incl_list$chain2 <- setNames(incl_chain2, paste0("b.", coef_names))
    fit <- structure(list(), class = "brmsfit")
    attr(fit, "sticky") <- list(
        inclusion_probs = incl_list,
        can_stick       = can_stick,
        unc_names       = unc_names
    )
    fit
}

# ──────────────────────────────────────────────────────────────────────────────
# inclusion_prob
# ──────────────────────────────────────────────────────────────────────────────

test_that("inclusion_prob.brmsfit returns named numeric in [0, 1] (single chain)", {
    fit <- make_mock_sticky_fit(c(0.9, 0.1, 0.8))
    ip  <- inclusion_prob(fit)
    expect_named(ip)
    expect_true(is.numeric(ip))
    expect_length(ip, 3L)
    expect_true(all(ip >= 0 & ip <= 1))
    expect_equal(names(ip), c("b.x1", "b.x2", "b.x3"))
})

test_that("inclusion_prob.brmsfit averages across chains", {
    fit <- make_mock_sticky_fit(c(0.8, 0.2), c(0.6, 0.4))
    ip  <- inclusion_prob(fit)
    expect_equal(ip, c(b.x1 = 0.7, b.x2 = 0.3))
})

test_that("inclusion_prob.brmsfit selects chain by index", {
    fit <- make_mock_sticky_fit(c(0.9, 0.1), c(0.5, 0.5))
    ip  <- inclusion_prob(fit, chain = 2L)
    expect_equal(ip, c(b.x1 = 0.5, b.x2 = 0.5))
})

test_that("inclusion_prob.brmsfit errors for missing chain", {
    fit <- make_mock_sticky_fit(c(0.9, 0.1))
    expect_error(inclusion_prob(fit, chain = 3L), "not found")
})

test_that("inclusion_prob errors on non-sticky brmsfit", {
    fit <- structure(list(), class = "brmsfit")
    expect_error(inclusion_prob(fit), "sticky")
})

test_that("inclusion_prob errors on non-brmsfit input", {
    expect_error(inclusion_prob(list()), class = "rlang_error")
})

# ──────────────────────────────────────────────────────────────────────────────
# median_probability_model
# ──────────────────────────────────────────────────────────────────────────────

test_that("median_probability_model returns names above threshold", {
    fit <- make_mock_sticky_fit(c(0.9, 0.3, 0.7))
    expect_warning(
        mpm <- median_probability_model(fit),
        "collinearity"
    )
    expect_equal(mpm, c("b.x1", "b.x3"))
})

test_that("median_probability_model with custom threshold", {
    fit <- make_mock_sticky_fit(c(0.9, 0.3, 0.7))
    mpm <- suppressWarnings(median_probability_model(fit, threshold = 0.8))
    expect_equal(mpm, "b.x1")
})

test_that("median_probability_model with threshold = 0 returns all names", {
    fit <- make_mock_sticky_fit(c(0.9, 0.01))
    mpm <- suppressWarnings(median_probability_model(fit, threshold = 0))
    expect_equal(mpm, c("b.x1", "b.x2"))
})

test_that("median_probability_model errors on invalid threshold", {
    fit <- make_mock_sticky_fit(c(0.5, 0.5))
    expect_error(median_probability_model(fit, threshold = 1.5), "threshold")
    expect_error(median_probability_model(fit, threshold = "a"), "threshold")
})

# ──────────────────────────────────────────────────────────────────────────────
# plot_inclusion
# ──────────────────────────────────────────────────────────────────────────────

test_that("plot_inclusion returns a ggplot object", {
    skip_if_not_installed("ggplot2")
    fit <- make_mock_sticky_fit(c(0.9, 0.1, 0.7))
    expect_warning(
        p <- plot_inclusion(fit),
        NA  # no warning expected from plot_inclusion itself
    )
    expect_s3_class(p, "gg")
})
