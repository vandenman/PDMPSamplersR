functions {
  // Pseudolikelihood numerator implemented without constructing full thresholds matrix.
//   real pseudolikelihood_numerator_v(
//       vector thresholds_0,                 // length = noThresholds
//       vector interactions_0,               // length = E
//       matrix suffstats,                    // P x P
//       array[] int seen,                    // size P
//       array[,] int threshold_counts_without_0, // P x K
//       int P,
//       int E,
//       array[] int edge_i,                        // length E
//       array[] int edge_j,                        // length E
//       array[] int threshold_start                // length P
//   ) {
//     real result = 0.0;
//     // thresholds contribution: sum counts * threshold
//     for (i in 1:P) {
//       int s = seen[i] - 1;
//       if (s > 0) {
//         for (u in 1:s) {
//           int tidx = threshold_start[i] + (u - 1);
//           result += threshold_counts_without_0[i, u] * thresholds_0[tidx];
//         }
//       }
//     }

//     // interactions contribution: sum_{i,j} interactions[i,j] * suffstats[i,j]
//     // original code used full symmetric interactions with zeros on diagonal.
//     // For interactions_0 storing only i<j edges, replicate the effect:
//     for (k in 1:E) {
//       int ii = edge_i[k];
//       int jj = edge_j[k];
//       // interactions[ii,jj] == interactions[jj,ii] == interactions_0[k]
//       // sum(interactions .* suffstats) originally sums both (ii,jj) and (jj,ii),
//       // so we add both entries explicitly
//       result += interactions_0[k] * (suffstats[ii, jj] + suffstats[jj, ii]);
//     }

//     return result;
//   }

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

  // Pseudolikelihood denominator: uses precomputed linpred (N x P) and
  // accesses thresholds via thresholds_0 + threshold_start mapping.
  // This is left as a single function so it can be parallelized later.
  real pseudolikelihood_denominator_v(
      vector exp_thresholds_0, // length noThresholds
      matrix linpred,      // N x P (X * interactions')
      array[] int seen,    // P
      int P,
      int N,
      array[] int threshold_start // length P
  ) {

    real result = 0.0;
    for (v in 1:N) {
      for (i in 1:P) {
        int s = seen[i] - 1;
        // vector[s+1] logits;
        if (s > 0) {
        //   logits[1] = 0.0;
          // fill logits[u] = thresholds[i,u] + u * linpred[v,i]
        //   for (u in 1:s) {
        //     int tidx = threshold_start[i] + (u - 1);
        //     logits[u+1] = thresholds_0[tidx] + u * linpred[v, i];
        //   }
        //   result += -log_sum_exp(logits);
          real temp0 = 1.0;
          real temp1 = exp(linpred[v, i]);
          real temp2 = temp1;
          for (u in 1:s) {
            int tidx = threshold_start[i] + (u - 1);
            temp0 += exp_thresholds_0[tidx] * temp1;
            temp1 *= temp2; // exp((u+1) * linpred[vi, i]) == exp(u * linpred[vi, i]) * exp(linpred[vi, i]);
          }
          result += -log(temp0);

        }// else {
          // seen[i] == 1 -> only 1 category => denominator contribution is -log(1) = 0
        //}
      }
    }
    return result;
  }
}

data {
  int N; // persons
  int P; // nodes
  int K; // maximum observed categories across all nodes
  matrix[N, P] X; // data (observations in {0,1,...})
  matrix[P, P] suffstats; // sufficient statistics
  array[P] int seen; // number of observed categories per node
  matrix[P, K] threshold_counts_without_0; // counts per node per threshold (u = 1..seen-1)
  real<lower = 0> prior_cauchy_scale;
  real<lower = 0> prior_threshold_alpha;
  real<lower = 0> prior_threshold_beta;
}

transformed data {
  // edge count and number of thresholds
  int E = P * (P - 1) / 2;
  int noThresholds = sum(seen) - P; // total (seen[i]-1) across i
  matrix[K, P] threshold_counts_without_0_t = threshold_counts_without_0'; // transpose for easier access

  array[P] int threshold_start;
  {
    int tidx = 1;
    for (i in 1:P) {
      threshold_start[i] = tidx;
      tidx += (seen[i] - 1);
    }
    // tidx-1 should equal noThresholds
  }
}

parameters {
  vector[noThresholds] thresholds_0; // all thresholds flattened (per-node)
  vector[E] interactions_0;          // unique i<j interactions flattened
}

transformed parameters {
  // Build symmetric interactions matrix once to use BLAS for linpred.
  // This avoids recomputing pairwise sums in the N*P loops.
  matrix[P, P] interactions;

  // initialize diagonal to zero
  for (i in 1:P)
    interactions[i, i] = 0.0;

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

  matrix[N, P] linpred = X * interactions;
  vector[noThresholds] exp_thresholds_0 = exp(thresholds_0);
}
// transformed parameters {
//   matrix[N, P] linpred = rep_matrix(0.0, N, P);
//   vector[noThresholds] exp_thresholds_0 = exp(thresholds_0);

//   // initialize to zero
//   {
//     int idx = 1;
//     for (i in 1:(P - 1)) {
//       for (j in (i + 1):P) {
//         // interactions_0[idx] = effect shared between i and j
//         real w = interactions_0[idx];
//         vector[N] Xi = X[, i];
//         vector[N] Xj = X[, j];

//         // update both columns at once using vectorized operations
//         linpred[, i] += w * Xj;
//         linpred[, j] += w * Xi;

//         idx += 1;
//       }
//     }
//   }
// }
model {
  // ---- priors (vectorized) ----
  target += sum(prior_threshold_alpha * thresholds_0 - (prior_threshold_alpha + prior_threshold_beta) * log1p_exp(thresholds_0));
  target += cauchy_lpdf(interactions_0 | 0, prior_cauchy_scale);

  // ---- pseudolikelihood numerator ----
  target += pseudolikelihood_numerator(
        thresholds_0,
        interactions,
        suffstats,
        threshold_counts_without_0_t
    );

    // target += pseudolikelihood_numerator2(
    //     thresholds_0,
    //     interactions_0,
    //     suffstats,
    //     threshold_counts_without_0_t,
    //     P
    // );

  // ---- pseudolikelihood denominator ----
  target += pseudolikelihood_denominator_v(
    exp_thresholds_0,
    linpred,
    seen,
    P,
    N,
    threshold_start
  );
}
