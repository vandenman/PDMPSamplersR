functions {
  int pdmp_get_subsample_size();
  int pdmp_get_subsample_index(int n);

  real omrf_person_lp(array[] int x, real interaction) {
    real eta_1 = interaction * x[2];
    real eta_2 = interaction * x[1];
    return x[1] * eta_1 - log1p_exp(eta_1) +
      x[2] * eta_2 - log1p_exp(eta_2);
  }
}
data {
  int<lower=1> N;
  array[N, 2] int<lower=0, upper=1> X;
  int<lower=0, upper=1> prior_only;
  real<lower=0> prior_sd;
}
parameters {
  real interaction;
}
model {
  interaction ~ normal(0, prior_sd);
  if (prior_only == 0) {
    int m = pdmp_get_subsample_size();
    if (m == 0) {
      for (n in 1:N) target += omrf_person_lp(X[n], interaction);
    } else {
      for (i in 1:m) {
        int n = pdmp_get_subsample_index(i);
        target += omrf_person_lp(X[n], interaction);
      }
    }
  }
}
