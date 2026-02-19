// ordinal_mrf_pseudolikelihood_glm_version.stan
data {
  int<lower=1> N; // persons
  int<lower=1> P; // nodes
  int<lower=1> K; // maximum possible categories (>=1), categories are 0..K-1 in usage
  matrix[N, P] X; // data (observations in {0,1,...})
  array[P] int seen; // number of observed categories per node (e.g. 1..K)
  real<lower = 0> prior_cauchy_scale;
  real<lower = 0> prior_threshold_alpha;
  real<lower = 0> prior_threshold_beta;
}

transformed data {
  // edge count and number of thresholds
  int E = P * (P - 1) / 2;
  int noThresholds = 0;
  for (i in 1:P) noThresholds += (seen[i] - 1);

  array[P] int threshold_start;
  {
    int tidx = 1;
    for (i in 1:P) {
      threshold_start[i] = tidx;
      tidx += (seen[i] - 1);
    }
    // tidx-1 should equal noThresholds
  }


  int num_cats = K + 1;
// category multipliers w = [0, 1, 2, ..., K]
  matrix[1, num_cats] w;
  w[1,1] = 0;
  for (m in 2:num_cats)
    w[1, m] = m - 1;

  array[N, P] int Xint;// = to_int(X + 1);
  for (n in 1:N) {
    for (p in 1:P) {
      Xint[n, p] = to_int(X[n, p] + 1);
    }
  }

}

parameters {
  vector[noThresholds] thresholds_0; // flattened per-node thresholds (these are mu_{i,u})
  vector[E] interactions_0;          // unique i<j interactions flattened
}

transformed parameters {
  // Build symmetric interactions matrix once
  matrix[P, P] interactions;
  for (i in 1:P)
    for (j in 1:P)
      interactions[i, j] = 0.0;

  {
    int idx = 1;
    for (i in 1:(P-1)) {
      for (j in (i+1):P) {
        interactions[i, j] = interactions_0[idx];
        interactions[j, i] = interactions_0[idx];
        idx += 1;
      }
    }
  }

  // compute linear predictor S = X * interactions (N x P)
  matrix[N, P] S = X * interactions;
}

model {
    // ---- priors (same as bgms) ----
    target += sum(prior_threshold_alpha * thresholds_0 -
      (prior_threshold_alpha + prior_threshold_beta) * log1p_exp(thresholds_0));
    target += cauchy_lpdf(interactions_0 | 0, prior_cauchy_scale);

    // ---- likelihood contribution vectorized over persons ----
    for (i in 1:P) {
        int s = seen[i] - 1; // number of thresholds for node i
        if (s > 0) {

            int tidx_start = threshold_start[i];   // 1-based index for thresholds_0

            // α_i must have length num_cats = K+1
            vector[num_cats] alpha_i;

            // set all to -∞ by default so unobserved categories get zero mass
            for (kk in 1:num_cats)
                alpha_i[kk] = negative_infinity();

            // category 0 has logit intercept 0
            alpha_i[1] = 0;

            // fill thresholds for categories 1..K
            for (u in 1:s)   // s = K
                alpha_i[u+1] = thresholds_0[tidx_start + u - 1];

            // Let me know if you find a way to avoid this (ugly) matrix construction
            // Build design matrix S_i (N×1)
            matrix[N,1] S_i;
            for (n in 1:N)
                S_i[n,1] = S[n,i];

            // GLM log-likelihood contribution for node i
            target += categorical_logit_glm_lpmf(Xint[, i] | S_i, alpha_i, w);

        }
    }
}
