functions {
  int pdmp_get_subsample_size();
  int pdmp_get_subsample_index(int n);

  real omrf_person_lp(array[] int x, row_vector x_real,
                      matrix thresholds, matrix interactions,
                      array[] int seen, int P) {
    real result = 0;
    row_vector[P] field = x_real * interactions;
    for (j in 1:P) {
      if (x[j] > 0) result += thresholds[j, x[j]];
      {
        real max_eta = 0;
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
        result += 2 * interactions[j, k] * x[j] * x[k];
    }
    return result;
  }
}

data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=1> K;
  int<lower=0, upper=1> prior_only;
  array[N, P] int<lower=0> X;
  array[P] int<lower=1> seen;
  real<lower=0> prior_interaction_sd;
  real<lower=0> prior_threshold_alpha;
  real<lower=0> prior_threshold_beta;
}

transformed data {
  int E = P * (P - 1) %/% 2;
  int noThresholds = sum(seen) - P;
  matrix[N, P] X_real;
  for (n in 1:N) for (j in 1:P) X_real[n, j] = X[n, j];
}

parameters {
  vector[noThresholds] thresholds_0;
  vector[E] interactions_0;
}

transformed parameters {
  matrix[P, K] thresholds = rep_matrix(negative_infinity(), P, K);
  matrix[P, P] interactions = rep_matrix(0, P, P);
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
  interactions_0 ~ normal(0, prior_interaction_sd);

  if (prior_only == 0) {
    int m = pdmp_get_subsample_size();
    if (m == 0) {
      for (n in 1:N)
        target += omrf_person_lp(X[n], X_real[n], thresholds,
                                 interactions, seen, P);
    } else {
      for (i in 1:m) {
        int n = pdmp_get_subsample_index(i);
        // Unscaled selected sum. Julia applies N / m exactly once.
        target += omrf_person_lp(X[n], X_real[n], thresholds,
                                 interactions, seen, P);
      }
    }
  }
}
