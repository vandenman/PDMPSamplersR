// ordinal_mrf_pseudolikelihood_glm_version.stan
// Variant of the user's model that uses categorical_logit_glm_lpmf per-node
// to replace the explicit denominator loop while keeping the paper's likelihood.
// (Based on the uploaded paper / your previous implementation.)

functions {
  real pseudolikelihood_numerator(
        vector thresholds_0,
        matrix interactions,
        matrix suffstats,
        matrix threshold_counts_without_0
    ) {

    real result = 0.0;

    result += dot_product(to_vector(threshold_counts_without_0), thresholds_0);
    result += dot_product(to_vector(interactions), to_vector(suffstats));

    return result;
  }

  // Keep this alternative numerator (edge-list form) if you want it
  real pseudolikelihood_numerator2(
        vector thresholds_0,
        vector interactions_0,
        matrix suffstats,
        matrix threshold_counts_without_0,
        int P) {

    real result = 0.0;

    result += dot_product(to_vector(threshold_counts_without_0), thresholds_0);
    // interactions contribution using edge list
    int k = 1;
    for (j in 1:P) {
      for (i in (j+1):P) {
        result += 2 * interactions_0[k] * suffstats[i, j];
        k = k + 1;
      }
    }
    return result;
  }
}

data {
  int<lower=1> N; // persons
  int<lower=1> P; // nodes
  int<lower=1> K; // maximum possible categories (>=1), categories are 0..K-1 in usage
  matrix[N, P] X; // data (observations in {0,1,...})
  matrix[P, P] suffstats; // sufficient statistics (same as before)
  array[P] int seen; // number of observed categories per node (e.g. 1..K)
  matrix[P, K] threshold_counts_without_0; // counts per node per threshold (u = 1..seen-1)
  real<lower = 0> prior_cauchy_scale;
  real<lower = 0> prior_threshold_alpha;
  real<lower = 0> prior_threshold_beta;
}

transformed data {
  // edge count and number of thresholds
  int E = P * (P - 1) / 2;
  int noThresholds = 0;
  for (i in 1:P) noThresholds += (seen[i] - 1); // total (seen[i]-1) across i

  matrix[K, P] threshold_counts_without_0_t = threshold_counts_without_0'; // transpose for ease

  array[P] int threshold_start;
  {
    int tidx = 1;
    for (i in 1:P) {
      threshold_start[i] = tidx;
      tidx += (seen[i] - 1);
    }
    // tidx-1 should equal noThresholds
  }

  // Build w = [0, 1, 2, ..., (K-1)] used as per-category multipliers in GLM
  // vector[K] w;
  // w[1] = 0;
  // for (m in 2:K) w[m] = m - 1;
  // matrix[K,1] w;
  // w[1,1] = 0;
  // for (k in 2:K) w[k,1] = k - 1;

  int num_cats = K + 1;
// category multipliers w = [0, 1, 2, ..., K]
  // matrix[num_cats, 1] w;
  // w[1,1] = 0;
  // for (m in 2:num_cats)
  //   w[m,1] = m - 1;
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
  vector[noThresholds] thresholds_0; // flattened per-node thresholds (these are mu_{i,u} in logit scale)
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
  matrix[N, P] S = X * interactions; // S[n,i] = sum_j X[n,j] * interactions[j,i]
  // S = X * interactions; // S[n,i] = sum_j X[n,j] * interactions[j,i]
  // matrix[P, N] St = interactions * X';
}

model {
  // ---- priors (same as your original) ----
  target += sum(prior_threshold_alpha * thresholds_0 - (prior_threshold_alpha + prior_threshold_beta) * log1p_exp(thresholds_0));
  target += cauchy_lpdf(interactions_0 | 0, prior_cauchy_scale);
  {
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

        // Build design matrix S_i (N×1)
        matrix[N,1] S_i;
        for (n in 1:N)
          S_i[n,1] = S[n,i];

        // GLM log-likelihood contribution for node i
        target += categorical_logit_glm_lpmf(Xint[, i] | S_i, alpha_i, w);

        // failed attempts to improve this further
        // target += categorical_logit_glm_lpmf(Xint[, i] | X, alpha_i, to_matrix( interactions[, i] )');
        // target += categorical_logit_glm_lpmf(Xint[, i] | to_matrix(S[, i]), alpha_i, w);
        // target += categorical_logit_glm_lpmf(Xint[, i] | St[i, ], alpha_i, w);
      }
    }
  }
}

// generated quantities {
//   int X_rep[N, P];
//   // generate posterior predictive (matching the pseudolikelihood form)
//   for (i in 1:P) {
//     int s = seen[i] - 1;
//     if (s > 0) {
//       int tidx_start = threshold_start[i];
//       vector[K] alpha_i;
//       for (kk in 1:K) alpha_i[kk] = negative_infinity();
//       alpha_i[1] = 0.0;
//       for (u in 1:s) alpha_i[u+1] = thresholds_0[tidx_start + u - 1];

//       for (n in 1:N) {
//         vector[K] logits;
//         // compute logits[1..s+1], fill others with very negative values
//         for (kk in 1:K) logits[kk] = negative_infinity();
//         for (kk in 1:(s+1)) logits[kk] = alpha_i[kk] + w[kk] * S[n, i]; // w[1]=0, w[2]=1,...
//         X_rep[n, i] = categorical_logit_rng(logits) - 1;
//       }
//     } else {
//       // only one category observed (category 0). Predict always 0
//       for (n in 1:N) X_rep[n, i] = 0;
//     }
//   }
// }
