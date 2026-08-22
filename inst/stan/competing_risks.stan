// Joint longitudinal model with competing risks (cause-specific Weibull hazards)
data {
  int<lower=1> N;
  int<lower=1> n_obs;
  int<lower=2> K;                        // number of causes
  vector[n_obs] y;
  vector[n_obs] time;
  array[n_obs] int<lower=1, upper=N> id;
  vector[N] surv_time;
  array[N] int<lower=0, upper=K> event; // 0 = censored, k = cause k
}
parameters {
  real beta_0;
  real beta_t;
  vector[2] alpha;                       // [value, slope] association
  matrix[K, 2] gamma_k;                  // cause-specific baseline params:
                                         // col 1 log-scale, col 2 log-shape
  real<lower=0> sigma_b0;
  real<lower=0> sigma_b1;
  real<lower=0> sigma_y;
  vector[N] z_b0;
  vector[N] z_b1;
}
model {
  vector[N] b0 = sigma_b0 * z_b0;
  vector[N] b1 = sigma_b1 * z_b1;

  // longitudinal
  {
    vector[n_obs] mu;
    for (j in 1:n_obs)
      mu[j] = beta_0 + b0[id[j]] + (beta_t + b1[id[j]]) * time[j];
    y ~ normal(mu, sigma_y);
  }
  z_b0 ~ std_normal();
  z_b1 ~ std_normal();

  for (i in 1:N) {
    real cv = beta_0 + b0[i] + (beta_t + b1[i]) * surv_time[i];
    real eta_v = alpha[1] * cv;
    real eta_s = alpha[2] * b1[i];
    for (k in 1:K) {
      real scale_k = exp(gamma_k[k, 1]);
      real shape_k = exp(gamma_k[k, 2]);
      real eta = eta_v + eta_s;
      if (event[i] == k)
        target += weibull_lpdf(surv_time[i] | shape_k,
                               scale_k * exp(-eta / shape_k));
      else if (event[i] == 0)
        target += weibull_lccdf(surv_time[i] | shape_k,
                                scale_k * exp(-eta / shape_k));
    }
  }

  // priors
  beta_0 ~ normal(0, 20);
  beta_t ~ normal(0, 5);
  alpha ~ normal(0, 1);
  to_vector(gamma_k) ~ normal(0, 3);
  sigma_b0 ~ exponential(1);
  sigma_b1 ~ exponential(1);
  sigma_y ~ exponential(1);
}
