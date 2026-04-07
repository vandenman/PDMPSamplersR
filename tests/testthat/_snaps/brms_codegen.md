# custom-family code snapshot: bernoulli fixed effects

    Code
      cat(result$ext_cpp)
    Output
      // generated with brms 2.23.0, modified by PDMPSamplersR 0.1.0
      functions {
          int pdmp_get_subsample_size();
        int pdmp_get_subsample_index(int n);
        matrix get_subsampled_Xc(matrix Xc_full);
        array[] int get_subsampled_int_array(array[] int arr);
          real subsampled_bernoulli_logit_lpmf(array[] int y, vector mu, int N_total) {
          int m = pdmp_get_subsample_size();
          real ll = 0;
          real scaling = 1.0 * N_total / m;
          for (i in 1:m) {
            int idx = pdmp_get_subsample_index(i);
            ll += bernoulli_logit_lpmf(y[idx] | mu[i]);
          }
          return ll * scaling;
        }
      }
      data {
        int<lower=1> N;  // total number of observations
        array[N] int Y;  // response variable
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int<lower=1> Kc;  // number of population-level effects after centering
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        vector[Kc] b;  // regression coefficients
        real Intercept;  // temporary intercept for centered predictors
      }
      transformed parameters {
        // prior contributions to the log posterior
        real lprior = 0;
        lprior += student_t_lpdf(Intercept | 3, 0, 2.5);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b;
          target += subsampled_bernoulli_logit_lpmf(Y | mu, N);
        }
        // priors including constants
        target += lprior;
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
      }

# custom-family code snapshot: poisson fixed effects

    Code
      cat(result$ext_cpp)
    Output
      // generated with brms 2.23.0, modified by PDMPSamplersR 0.1.0
      functions {
          int pdmp_get_subsample_size();
        int pdmp_get_subsample_index(int n);
        matrix get_subsampled_Xc(matrix Xc_full);
        array[] int get_subsampled_int_array(array[] int arr);
          real subsampled_poisson_log_lpmf(array[] int y, vector mu, int N_total) {
          int m = pdmp_get_subsample_size();
          real ll = 0;
          real scaling = 1.0 * N_total / m;
          for (i in 1:m) {
            int idx = pdmp_get_subsample_index(i);
            ll += poisson_log_lpmf(y[idx] | mu[i]);
          }
          return ll * scaling;
        }
      }
      data {
        int<lower=1> N;  // total number of observations
        array[N] int Y;  // response variable
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int<lower=1> Kc;  // number of population-level effects after centering
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        vector[Kc] b;  // regression coefficients
        real Intercept;  // temporary intercept for centered predictors
      }
      transformed parameters {
        // prior contributions to the log posterior
        real lprior = 0;
        lprior += student_t_lpdf(Intercept | 3, 0.7, 2.5);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b;
          target += subsampled_poisson_log_lpmf(Y | mu, N);
        }
        // priors including constants
        target += lprior;
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
      }

# custom-family code snapshot: gaussian fixed effects

    Code
      cat(result$ext_cpp)
    Output
      // generated with brms 2.23.0, modified by PDMPSamplersR 0.1.0
      functions {
          int pdmp_get_subsample_size();
        int pdmp_get_subsample_index(int n);
        matrix get_subsampled_Xc(matrix Xc_full);
        array[] int get_subsampled_int_array(array[] int arr);
          real subsampled_gaussian_lpdf(vector y, vector mu, real sigma, int N_total) {
          int m = pdmp_get_subsample_size();
          real ll = 0;
          real scaling = 1.0 * N_total / m;
          for (i in 1:m) {
            int idx = pdmp_get_subsample_index(i);
            ll += normal_lpdf(y[idx] | mu[i], sigma);
          }
          return ll * scaling;
        }
      }
      data {
        int<lower=1> N;  // total number of observations
        vector[N] Y;  // response variable
        int<lower=1> K;  // number of population-level effects
        matrix[N, K] X;  // population-level design matrix
        int<lower=1> Kc;  // number of population-level effects after centering
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
        matrix[N, Kc] Xc;  // centered version of X without an intercept
        vector[Kc] means_X;  // column means of X before centering
        for (i in 2:K) {
          means_X[i - 1] = mean(X[, i]);
          Xc[, i - 1] = X[, i] - means_X[i - 1];
        }
      }
      parameters {
        vector[Kc] b;  // regression coefficients
        real Intercept;  // temporary intercept for centered predictors
        real<lower=0> sigma;  // dispersion parameter
      }
      transformed parameters {
        // prior contributions to the log posterior
        real lprior = 0;
        lprior += student_t_lpdf(Intercept | 3, -0.1, 2.5);
        lprior += student_t_lpdf(sigma | 3, 0, 2.5)
          - 1 * student_t_lccdf(0 | 3, 0, 2.5);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b;
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
      }
      generated quantities {
        // actual population-level intercept
        real b_Intercept = Intercept - dot_product(means_X, b);
      }

# rewrite snapshot: rewrite_mu_init on FE model

    Code
      cat(result)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + Xc * b;
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;

# rewrite snapshot: rewrite_fe_accumulation on FE model

    Code
      cat(result)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + get_subsampled_Xc(Xc) * b;
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;

# rewrite snapshot: rewrite_spline_matrices on spline model

    Code
      cat(result)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b + get_subsampled_Xc(Xs) * bs + get_subsampled_Xc(Zs_1_1) * s_1_1;
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(zs_1_1);

# rewrite snapshot: rewrite_gp_indexing on GP model

    Code
      cat(result)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          vector[Nsubgp_1] gp_pred_1 = gp_exp_quad(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + gp_pred_1[get_subsampled_int_array(Jgp_1)];
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(zgp_1);

# rewrite snapshot: rewrite_re_loops on RE intercept model

    Code
      cat(result)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b;
          for (i in 1:m_sub) {
            int n = pdmp_get_subsample_index(i);
            // add more terms to the linear predictor
            mu[i] += r_1_1[J_1[n]] * Z_1_1[n];
          }
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);

# rewrite snapshot: rewrite_re_loops on RE slope model

    Code
      cat(result)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b;
          for (i in 1:m_sub) {
            int n = pdmp_get_subsample_index(i);
            // add more terms to the linear predictor
            mu[i] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
          }
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(to_vector(z_1));

# rewrite snapshot: rewrite_re_loops on multi-group RE model

    Code
      cat(result)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] mu = rep_vector(0.0, N);
          mu += Intercept + Xc * b;
          for (i in 1:m_sub) {
            int n = pdmp_get_subsample_index(i);
            // add more terms to the linear predictor
            mu[i] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
          }
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);
        target += std_normal_lpdf(z_2[1]);

# model block snapshot: gaussian y ~ x

    Code
      cat(block)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b;
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;

# model block snapshot: gaussian y ~ x + s(z)

    Code
      cat(block)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b + get_subsampled_Xc(Xs) * bs + get_subsampled_Xc(Zs_1_1) * s_1_1;
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(zs_1_1);

# model block snapshot: gaussian y ~ gp(x)

    Code
      cat(block)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          vector[Nsubgp_1] gp_pred_1 = gp_exp_quad(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + gp_pred_1[get_subsampled_int_array(Jgp_1)];
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(zgp_1);

# model block snapshot: gaussian y ~ x + (1 | group)

    Code
      cat(block)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b;
          for (i in 1:m_sub) {
            int n = pdmp_get_subsample_index(i);
            // add more terms to the linear predictor
            mu[i] += r_1_1[J_1[n]] * Z_1_1[n];
          }
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);

# model block snapshot: bernoulli y ~ x + (1 | group)

    Code
      cat(block)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b;
          for (i in 1:m_sub) {
            int n = pdmp_get_subsample_index(i);
            // add more terms to the linear predictor
            mu[i] += r_1_1[J_1[n]] * Z_1_1[n];
          }
          target += subsampled_bernoulli_logit_lpmf(Y | mu, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);

# model block snapshot: gaussian y ~ x + (1 + x | group)

    Code
      cat(block)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b;
          for (i in 1:m_sub) {
            int n = pdmp_get_subsample_index(i);
            // add more terms to the linear predictor
            mu[i] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
          }
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(to_vector(z_1));

# model block snapshot: gaussian y ~ x + (1 | g1) + (1 | g2)

    Code
      cat(block)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b;
          for (i in 1:m_sub) {
            int n = pdmp_get_subsample_index(i);
            // add more terms to the linear predictor
            mu[i] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
          }
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);
        target += std_normal_lpdf(z_2[1]);

# model block snapshot: gaussian y ~ x + s(z) + (1 | group)

    Code
      cat(block)
    Output
      
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          int m_sub = pdmp_get_subsample_size();
          vector[m_sub] mu = rep_vector(0.0, m_sub);
          mu += Intercept + get_subsampled_Xc(Xc) * b + get_subsampled_Xc(Xs) * bs + get_subsampled_Xc(Zs_1_1) * s_1_1;
          for (i in 1:m_sub) {
            int n = pdmp_get_subsample_index(i);
            // add more terms to the linear predictor
            mu[i] += r_1_1[J_1[n]] * Z_1_1[n];
          }
          target += subsampled_gaussian_lpdf(Y | mu, sigma, N);
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(zs_1_1);
        target += std_normal_lpdf(z_1[1]);

