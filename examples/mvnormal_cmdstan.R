# renv::install(".", prompt = FALSE)
library(PDMPSamplersR)
library(tibble)
library(ggplot2)
library(patchwork)
library(fs)

theme_set(
  theme_bw(base_size = 32) +
    theme(geom = element_geom(pointshape = 21, pointsize = 15, fill = scales::alpha("grey", .5)))
)

stan_model_dir <- path("stan", "models")
stan_data_dir  <- path("stan", "data")

set.seed(123)
d <- 5
mean_vec <- rnorm(d, 0, 2)
cov_matrix <- stats::rWishart(1, df = d + 2, Sigma = diag(d))[, , 1]
stan_data <- list(N = 5, mu = mean_vec, sigma = cov_matrix)

data_path <- path(stan_data_dir, "mvnormal_data.json")
write_stan_json(stan_data, data_path)

model_path <- path(stan_model_dir, "mvnormal.stan")
result <- PDMPSamplersR::pdmp_sample_from_stanmodel(model_path, data_path, flow = "ZigZag", T = 10000, flow_mean = mean_vec, flow_cov = cov_matrix)

result
cov(result[[1]]) - cov_matrix
colMeans(result[[1]]) - mean_vec

est_mean <- colMeans(result[[1]])
est_cov  <- cov(result[[1]])

tib <- tibble(
  truth       = c(mean_vec, cov_matrix[lower.tri(cov_matrix, diag = TRUE)]),
  estimate_zz = c(est_mean, est_cov[lower.tri(est_cov, diag = TRUE)]),
  parameter   = c(paste0("mu[", 1:d, "]"), unlist(lapply(1:d, function(i) paste0("sigma[", i:d, ", ", i, "]")))),
  what        = rep(c("mu", "sigma"), times = c(d, d * (d + 1) / 2))
)

ggplot(tib, aes(x = truth, y = estimate_zz)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .8) +
  geom_point(shape = 21, size = 12, fill = "grey", alpha = .65) +
  facet_wrap(~what, scales = "free") +
  theme_bw(base_size = 30) +
  labs(title = "ZigZag PDMP", x = "True value", y = "Estimated value")

result2 <- PDMPSamplersR::pdmp_sample_from_stanmodel(model_path, data_path, flow = "ZigZag", alg = "GridThinningStrategy", T = 10000)
result2
cov(result2[[1]]) - cov_matrix
colMeans(result2[[1]]) - mean_vec

est_mean <- colMeans(result2[[1]])
est_cov  <- cov(result2[[1]])

tib <- tibble(
  truth       = c(mean_vec, cov_matrix[lower.tri(cov_matrix, diag = TRUE)]),
  estimate_zz = c(est_mean, est_cov[lower.tri(est_cov, diag = TRUE)]),
  parameter   = c(paste0("mu[", 1:d, "]"), unlist(lapply(1:d, function(i) paste0("sigma[", i:d, ", ", i, "]")))),
  what        = rep(c("mu", "sigma"), times = c(d, d * (d + 1) / 2))
)

ggplot(tib, aes(x = truth, y = estimate_zz)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .8) +
  geom_point(shape = 21, size = 12, fill = "grey", alpha = .65) +
  facet_wrap(~what, scales = "free") +
  theme_bw(base_size = 30) +
  labs(title = "ZigZag PDMP", x = "True value", y = "Estimated value")


set.seed(123)
d <- 12
mean_vec <- rep(0, d)
cov_matrix <- stats::rWishart(1, df = d + 2, Sigma = diag(d))[, , 1]
cov_matrix <- diag(diag(cov_matrix))
stan_data <- list(N = d, mu = mean_vec, sigma = cov_matrix)

data_path <- path(stan_data_dir, "mvnormal_spike_and_slab_data.json")
write_stan_json(stan_data, data_path)


prior_incl_prob <- seq(0.1, 0.9, length.out = d)
result3 <- PDMPSamplersR::pdmp_sample_from_stanmodel(
  model_path, data_path,
  flow = "ZigZag", T = 10000,
  alg = "GridThinningStrategy",
  sticky = TRUE, can_stick = rep(TRUE, d),
  model_prior = bernoulli(prior_incl_prob),
  parameter_prior = dnorm(rep(0, d), mean_vec, sqrt(diag(cov_matrix)))
)

tib <- tibble(true = prior_incl_prob, est  = colMeans(result3[[1]] != 0))
ggplot(tib, aes(x = true, y = est)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .8) +
  geom_point() +
  labs(title = "ZigZag PDMP with variable selection", x = "True inclusion probability", y = "Estimated inclusion probability")



# BIASED! needs a function because the conditionals change!
set.seed(123)
d <- 12
mean_vec <- rep(0, d)
cov_matrix <- stats::rWishart(1, df = d + 2, Sigma = diag(d))[, , 1]
stan_data <- list(N = d, mu = mean_vec, sigma = cov_matrix)

prior_incl_prob <- seq(0.1, 0.9, length.out = d)

data_path <- path(stan_data_dir, "mvnormal_spike_and_slab_data2.json")
write_stan_json(stan_data, data_path)

result4 <- PDMPSamplersR::pdmp_sample_from_stanmodel(
  model_path, data_path,
  flow = "ZigZag", T = 10000,
  alg = "GridThinningStrategy",
  sticky = TRUE, can_stick = rep(TRUE, d),
  model_prior = bernoulli(prior_incl_prob),
  parameter_prior = dnorm(rep(0, d), mean_vec, sqrt(diag(cov_matrix)))
)

tib <- tibble(true = prior_incl_prob, est  = colMeans(result4[[1]] != 0))
ggplot(tib, aes(x = true, y = est)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .8) +
  geom_point() +
  labs(title = "ZigZag PDMP with variable selection", x = "True inclusion probability", y = "Estimated inclusion probability")

# results are biased because covariance matrix is not diagonal

# different model prior: beta-binomial
set.seed(123)
d <- 7
mean_vec <- rep(0, d)
cov_matrix <- stats::rWishart(1, df = d + 2, Sigma = diag(d))[, , 1]
cov_matrix <- diag(diag(cov_matrix))
stan_data <- list(N = d, mu = mean_vec, sigma = cov_matrix)

model_prior = betabernoulli(1, 1)
prior_incl_prob <- extraDistr::dbbinom(0:d, size = 1, alpha = model_prior$a, beta = model_prior$b)
beta_binomial_pmf <- function(k, d, a, b, log = FALSE, individual_model = FALSE) {
  # Inputs:
  # - k: number of ones (successes)
  # - d: total number of trials (length of binary vector)
  # - a, b: Beta prior parameters

  # Log binomial coefficient: log(C(d, k)) = lfactorial(d) - lfactorial(k) - lfactorial(d - k)
  log_binomial_coeff <- lfactorial(d) - lfactorial(k) - lfactorial(d - k)

  # Log Beta function terms: log(B(a+k, b+d-k)) and log(B(a, b))
  log_beta_numerator <- lbeta(a + k, b + d - k)
  log_beta_denominator <- lbeta(a, b)

  # Combine in log space and exponentiate
  log_prob <- (!individual_model) * log_binomial_coeff + log_beta_numerator - log_beta_denominator
  if (log) {
     return(log_prob)
  }
  exp(log_prob)
}
# beta_binomial_pmf(0:d, d, model_prior$a, model_prior$b)
prior_incl_prob <- rep(model_prior$a / (model_prior$a + model_prior$b), d)

data_path <- path(stan_data_dir, "mvnormal_spike_and_slab_data_beta_binomial.json")
write_stan_json(stan_data, data_path)

result5 <- PDMPSamplersR::pdmp_sample_from_stanmodel(
  model_path, data_path,
  flow = "ZigZag", T = 10000,
  alg = "GridThinningStrategy",
  sticky = TRUE, can_stick = rep(TRUE, d),
  model_prior = bernoulli(prior_incl_prob),
  parameter_prior = dnorm(rep(0, d), mean_vec, sqrt(diag(cov_matrix)))
)

tib_incl <- tibble(true = prior_incl_prob, est  = colMeans(result5[[1]] != 0))
plt_incl <- ggplot(tib_incl, aes(x = true, y = est)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .8) +
  geom_point() +
  labs(title = "Inclusion probabilities", x = "True inclusion probability", y = "Estimated inclusion probability")

model_frequencies <- table(apply(1L * (result5[[1]] != 0), 1L, paste0, collapse = ""))
model_probabilities <- model_frequencies / sum(model_frequencies)

# the no. 1s in the models
model_sizes <- nchar(gsub("0", "", names(model_frequencies), fixed = TRUE))

all_model_probabilities <- beta_binomial_pmf(0:d, d, model_prior$a, model_prior$b, log = TRUE, individual_model = TRUE)
expected_model_probabilities <- all_model_probabilities[model_sizes + 1]

tib_model <- tibble(
  true_log_model_probs = expected_model_probabilities,
  est_log_model_probs  = log(unname(c(model_probabilities))),
)
plt_model <- ggplot(tib_model, aes(x = true_log_model_probs, y = est_log_model_probs)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .8) +
  geom_point() +
  labs(title = "Model probabilities", x = "True model probability", y = "Estimated model probability")

plt_incl + plt_model



set.seed(123)
n <- 200
d <- 10
d_on <- floor(0.2 * d)
idx_on <- sample(1:d, d_on, FALSE)
intercept <- -3
beta <- rep(0, d)
beta[idx_on] <- rnorm(d_on, 0, 2)
X <- matrix(rnorm(n * d), n, d)
X <- scale(X)
linpred <- intercept + X %*% beta
y <- c(runif(n, 0, 1) <= plogis(linpred))
fit <- stats::glm(y ~ X, family = binomial())
coef(fit)
cor(coef(fit), c(intercept, beta))
sqrt(mean((coef(fit) - c(intercept, beta))^2))

stan_data <- list(N = n, D = d, X = X, y = y, sd_prior = .1)
data_path <- path(stan_data_dir, "logistic_regression_data.json")
model_path <- path(stan_model_dir, "logistic_regression.stan")
write_stan_json(stan_data, data_path)

result_lr <- PDMPSamplersR::pdmp_sample_from_stanmodel(
  model_path, data_path,
  flow = "ZigZag", T = 10000,
  alg = "GridThinningStrategy",
  sticky = TRUE,
  can_stick = c(FALSE, rep(TRUE, d)), # first parameter is intercept, which is always included
  model_prior = bernoulli(0.5),
  parameter_prior = dnorm(rep(0, d+1))
)

tib_lr <- tibble(
  true   = rep(c(intercept, beta), 2),
  est    = c(colMeans(result_lr[[1]]), coef(fit)),
  how    = rep(c("Bayes", "Frequentist"), each = d + 1),
  p_incl = c(colMeans(result_lr[[1]] != 0), rep(1, d+1)),
  param  = rep(c("intercept", rep("beta", d)), 2)
)
plt_recovery <- ggplot(tib_lr, aes(x = true, y = est, fill = how, shape = param)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .8) +
  geom_point(alpha = .5) +
  scale_shape_manual(values = c(21, 24)) +
  facet_wrap(~how) +
  labs(title = "Parameter recovery", x = "True parameter value", y = "Estimated parameter value")

plt_incl_vs_value <- subset(tib_lr, how == "Bayes" & param == "beta") |>
  ggplot(aes(x = est, y = p_incl)) +
  geom_point() +
  labs(title = "Inclusion probability vs. effect size", x = "Estimated posterior mean", y = "Inclusion probability")

plt_recovery / plt_incl_vs_value

# fs::file_delete("stan/models/logistic_regression_model.so")
sd_prior <- 10
stan_data <- list(N = n, D = d, X = X, y = y, sd_prior = sd_prior)
data_path <- path(stan_data_dir, "logistic_regression_data_sd_10.json")
write_stan_json(stan_data, data_path)

result_lr_sd10 <- PDMPSamplersR::pdmp_sample_from_stanmodel(
  model_path, data_path,
  flow = "ZigZag", T = 20000,
  alg = "GridThinningStrategy",
  sticky = TRUE,
  can_stick = c(FALSE, rep(TRUE, d)), # first parameter is intercept, which is always included
  model_prior = bernoulli(0.5),
  parameter_prior = dnorm(rep(0, d+1), 0, sd_prior)
)

# lib = "stan/models/logistic_regression_model.so"
# JuliaCall::julia_eval(sprintf("begin
# import Libdl
# Libdl.dlclose(%s)
# end", lib))

tib_lr <- tibble(
  true   = rep(c(intercept, beta), 2),
  est    = c(colMeans(result_lr_sd10[[1]]), coef(fit)),
  how    = rep(c("Bayes", "Frequentist"), each = d + 1),
  p_incl = c(colMeans(result_lr_sd10[[1]] != 0), rep(1, d+1)),
  param  = rep(c("intercept", rep("beta", d)), 2)
)
plt_recovery <- ggplot(tib_lr, aes(x = true, y = est, fill = how, shape = param)) +
  geom_abline(slope = 1, intercept = 0, color = "grey", alpha = .8) +
  geom_point(alpha = .5) +
  scale_shape_manual(values = c(21, 24)) +
  facet_wrap(~how) +
  labs(title = "Parameter recovery", x = "True parameter value", y = "Estimated parameter value")

plt_incl_vs_value <- subset(tib_lr, how == "Bayes" & param == "beta") |>
  ggplot(aes(x = est, y = p_incl)) +
  geom_point() +
  labs(title = "Inclusion probability vs. effect size", x = "Estimated posterior mean", y = "Inclusion probability")

plt_recovery / plt_incl_vs_value


# t-test
library(PDMPSamplersR)
set.seed(123)
nx <- 30
ny <- 35
mux <- 0
muy <- 0.5
sigma2 <- 2.1

x1 <- rnorm(nx, mux, sqrt(sigma2))
x2 <- rnorm(ny, muy, sqrt(sigma2))

rscale <- 1.0
analytic_results <- BayesFactor::ttestBF(x = x1, y = x2, rscale = rscale)

stan_data <- list(nx = nx, ny = ny, x = x1, y = x2, rscale = rscale)
data_path <- path(stan_data_dir, "ttest_data.json")
model_path <- path(stan_model_dir, "ttest.stan")
write_stan_json(stan_data, data_path)

# ignore this for now
result_ttest <- PDMPSamplersR::pdmp_sample_from_stanmodel(
  model_path, data_path,
  flow = "ZigZag", T = 100000,
  alg = "GridThinningStrategy",
  grid_n = 30,
  sticky = TRUE,
  can_stick = c(FALSE, TRUE, FALSE), # grand mean, delta, sigma
  model_prior = bernoulli(0.5),
  parameter_prior = rep(dcauchy(0, 0, rscale), 3)
)

prior_inclusion_prob <- 0.5
posterior_inclusion_prob <- mean(result_ttest$samples[, 2] != 0)
prior_odds <- prior_inclusion_prob / (1 - prior_inclusion_prob)
posterior_odds <- posterior_inclusion_prob / (1 - posterior_inclusion_prob)
bayes_factor_10 <- posterior_odds / prior_odds

# TODO: close but not close enough!
c(
  pdmp_bf = bayes_factor_10,
  analytic_bf = BayesFactor::extractBF(analytic_results, onlybf = TRUE)
)

# TODO: should be automatic with bridgestan!
# mean(exp(result_ttest$samples[, 3]))
