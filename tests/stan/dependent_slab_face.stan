data {
  int<lower=0, upper=1> face_mode;
  array[2] int<lower=0, upper=1> active;
  real y;
}
transformed data {
  vector[2] slab_mean = [0.2, -0.1]';
  matrix[2, 2] slab_cov = [[1.1, 0.35], [0.35, 0.8]];
}
parameters {
  vector[2] beta;
  real nuisance;
}
model {
  nuisance ~ normal(1, 0.7);
  y ~ normal(beta[1] - 0.5 * beta[2] + nuisance, 0.8);
  if (face_mode == 0) {
    beta ~ multi_normal(slab_mean, slab_cov);
  } else if (active[1] == 1 && active[2] == 1) {
    beta ~ multi_normal(slab_mean, slab_cov);
  } else if (active[1] == 1) {
    beta[1] ~ normal(slab_mean[1], sqrt(slab_cov[1, 1]));
  } else if (active[2] == 1) {
    beta[2] ~ normal(slab_mean[2], sqrt(slab_cov[2, 2]));
  }
}
