functions {
  int pdmp_get_subsample_size();
  int pdmp_get_subsample_index(int n);
}
data {
  int<lower=1> N;
  vector[N] y;
  int<lower=0, upper=1> prior_only;
}
parameters {
  real beta;
}
model {
  beta ~ normal(0, 2);
  if (prior_only == 0) {
    int m = pdmp_get_subsample_size();
    if (m == 0) {
      y ~ normal(beta, 1);
    } else {
      for (i in 1:m) {
        int n = pdmp_get_subsample_index(i);
        // Deliberately unscaled: Julia applies N / m.
        target += normal_lpdf(y[n] | beta, 1);
      }
    }
  }
}
