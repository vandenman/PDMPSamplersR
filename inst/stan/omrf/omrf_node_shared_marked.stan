functions {
  int pdmp_get_subsample_size();
  int pdmp_get_subsample_index(int n);

  real omrf_person_lp(array[] int x, row_vector x_real,
                      matrix thresholds, matrix interactions,
                      array[] int seen, int P) {
    real result = 0.0;
    row_vector[P] field = x_real * interactions;
    for (j in 1:P) {
      if (x[j] > 0) result += thresholds[j, x[j]];
      {
        real max_eta = 0.0;
        real denominator;
        if (seen[j] > 1) {
          for (u in 1:(seen[j] - 1))
            max_eta = fmax(max_eta, thresholds[j, u] + u * field[j]);
        }
        denominator = exp(-max_eta);
        if (seen[j] > 1) {
          for (u in 1:(seen[j] - 1))
            denominator += exp(thresholds[j, u] + u * field[j] - max_eta);
        }
        result -= max_eta + log(denominator);
      }
    }
    if (P > 1) {
      for (j in 1:(P - 1)) for (k in (j + 1):P)
        result += 2.0 * interactions[j, k] * x[j] * x[k];
    }
    return result;
  }
}

data {
  int<lower=1> N;
  int<lower=2> P;
  int<lower=2> K;
  int<lower=0, upper=1> prior_only;
  array[N, P] int<lower=0> X;
  array[P] int<lower=1> seen;
  real<lower=0> prior_threshold_alpha;
  real<lower=0> prior_threshold_beta;
  real<lower=0> slab_base_scale;
  real<lower=0> prior_log_tau_sd;
  real<lower=0> prior_log_lambda_sd;
}

transformed data {
  int E = (P * (P - 1)) %/% 2;
  int noThresholds = sum(seen) - P;
  matrix[N, P] X_real;
  for (n in 1:N) for (j in 1:P) X_real[n, j] = X[n, j];
}

parameters {
  vector[noThresholds] thresholds_0;
  vector[E] interactions_0;
  real log_tau;
  vector[P] log_lambda;
}

transformed parameters {
  matrix[P, K] thresholds = rep_matrix(negative_infinity(), P, K);
  matrix[P, P] interactions = rep_matrix(0.0, P, P);
  {
    int pos = 1;
    for (j in 1:P) {
      if (seen[j] > 1) for (u in 1:(seen[j] - 1)) {
        thresholds[j, u] = thresholds_0[pos];
        pos += 1;
      }
    }
    pos = 1;
    if (P > 1) for (j in 1:(P - 1)) for (k in (j + 1):P) {
      interactions[j, k] = interactions_0[pos];
      interactions[k, j] = interactions_0[pos];
      pos += 1;
    }
  }
}

model {
  target += sum(prior_threshold_alpha * thresholds_0 -
                (prior_threshold_alpha + prior_threshold_beta) *
                log1p_exp(thresholds_0));
  log_tau ~ normal(0, prior_log_tau_sd);
  log_lambda ~ normal(0, prior_log_lambda_sd);
  {
    int edge = 1;
    for (j in 1:(P - 1)) for (k in (j + 1):P) {
      real log_sd = log(slab_base_scale) + log_tau
        + 0.5 * log_lambda[j] + 0.5 * log_lambda[k];
      interactions_0[edge] ~ normal(0, exp(log_sd));
      edge += 1;
    }
  }

  if (prior_only == 0) {
    int m = pdmp_get_subsample_size();
    if (m == 0) {
      for (n in 1:N)
        target += omrf_person_lp(X[n], X_real[n], thresholds,
                                 interactions, seen, P);
    } else {
      for (i in 1:m) {
        int n = pdmp_get_subsample_index(i);
        target += omrf_person_lp(X[n], X_real[n], thresholds,
                                 interactions, seen, P);
      }
    }
  }
}
