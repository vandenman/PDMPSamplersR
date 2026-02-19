functions{
  real pseudolikelihood_numerator(
        matrix thresholds,
        matrix interactions,
        matrix suffstats,
        array[] int seen,
        array[,] int threshold_counts_without_0,
        int P) {

    real result = 0.0;
    for (i in 1:P) {
      for (u in 1:(seen[i] - 1)) {
        result += threshold_counts_without_0[i, u] * thresholds[i, u];
      }
    }

    result += sum(interactions .* suffstats);

    return result;
  }

  real pseudolikelihood_denominator(
      matrix thresholds,
      matrix interactions,
      matrix suffstats,
      array[] int seen,
      array[,] int X,
      int P, int N
  ) {
    real result = 0.0;
    for (v in 1:N) {
      for (i in 1:P) {
        real temp1 = 0.0;
        for (j in 1:P) {
          // note that interactions[i, i] = 0
          temp1 += interactions[i, j] * X[v, j];
        }
        real temp2 = 1.0;
        for (u in 1:(seen[i]-1)) {
          temp2 += exp(thresholds[i, u] + u * temp1);
        }
        result += -log(temp2);
      }
    }
    return result;
  }

}
data {
  int N; // persons
  int P; // nodes
  int K; // maximum observed categories across all nodes
  array[N, P] int X; // data
  matrix[P, P] suffstats; // sufficient statistics
  array[P] int seen; // no observed categories per node
  array[P, K] int threshold_counts_without_0; // no.
  real<lower = 0> prior_cauchy_scale;
  real<lower = 0> prior_threshold_alpha;
  real<lower = 0> prior_threshold_beta;
}
transformed data {
  int E = P * (P - 1) / 2; // edges
  int noThresholds = sum(seen) - P;
}
parameters {
  vector[noThresholds] thresholds_0;
  vector[E] interactions_0;
}
transformed parameters {
  // this is not efficient, just easy to write up
  matrix[P, K] thresholds;
  matrix[P, P] interactions;
  // define new local scope for idx
  {
    int idx = 1;
    for (i in 1:P) {
      for (j in 1:(seen[i] - 1)) {
        thresholds[i, j] = thresholds_0[idx];
        idx += 1;
      }
      for (j in (seen[i]) : K) {
        thresholds[i, j] = negative_infinity();
      }
    }
    idx = 1;
    for (i in 1:(P-1)) {
      interactions[i, i] = 0;
      for (j in (i+1):P) {
        interactions[i, j] = interactions_0[idx];
        interactions[j, i] = interactions_0[idx];
        idx += 1;
      }
    }
    interactions[P, P] = 0;
  }
}
model {
  // priors
  // for (i in 1:noThresholds) {
  //   target += prior_threshold_alpha * thresholds_0[i] -
  //     (prior_threshold_alpha + prior_threshold_beta) * log1p_exp(thresholds_0[i]);
  // }
  target += sum(prior_threshold_alpha * thresholds_0 - (prior_threshold_alpha + prior_threshold_beta) * log1p_exp(thresholds_0));
  target += cauchy_lpdf(interactions_0 | 0, prior_cauchy_scale);

  // pseudolikelihood numerator
  for (i in 1:P) {
    for (u in 1:(seen[i] - 1)) {
      target += threshold_counts_without_0[i, u] * thresholds[i, u];
    }
  }
  target += sum(interactions .* suffstats);
  //pseudolikelihood denominator
  // vector[K] temp;
  for (v in 1:N) {
    for (i in 1:P) {
      real temp1 = 0.0;
      for (j in 1:P) {
        // note that interactions[i, i] = 0
        temp1 += interactions[i, j] * X[v, j];
      }
      real temp2 = 1.0;
      for (u in 1:(seen[i]-1)) {
        temp2 += exp(thresholds[i, u] + u * temp1);
      }
      target += -log(temp2);
    }
  }
}

