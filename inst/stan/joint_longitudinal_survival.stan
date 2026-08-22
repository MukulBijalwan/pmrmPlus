// Joint longitudinal-survival model (Weibull hazard, current-value association)
data {
  int<lower=1> N;                       // subjects
  int<lower=1> n_obs;                   // total observations
  int<lower=1> P;                       // fixed effects
  vector[n_obs] y;
  vector[n_obs] time;
  array[n_obs] int<lower=1, upper=N> id;
  vector[N] surv_time;
  array[N] int<lower=0, upper=1> event; // 1 = event, 0 = censored
  vector[N] trt;
}
parameters {
  vector[P] beta;
  array[2] vector[N] b_std;             // non-centered random intercept/slope
  real<lower=0> sigma_b0;
  real<lower=0> sigma_b1;
  real<lower=0> sigma_y;
  real alpha_trt;
  real alpha_assoc;
  real<lower=0> lambda_weibull;
  real<lower=0> nu_weibull;
}
transformed parameters {
  vector[n_obs] yhat;
  for (j in 1:n_obs) {
    yhat[j] = beta[1] + sigma_b0 * b_std[1][id[j]] +
              (beta[2] + sigma_b1 * b_std[2][id[j]]) * time[j];
  }
}
model {
  // longitudinal likelihood
  y ~ normal(yhat, sigma_y);

  // non-centered random effects
  b_std[1] ~ std_normal();
  b_std[2] ~ std_normal();

  // survival: Weibull with current-value association
  for (i in 1:N) {
    real cv = beta[1] + sigma_b0 * b_std[1][i] +
              (beta[2] + sigma_b1 * b_std[2][i]) * surv_time[i];
    real eta = alpha_trt * trt[i] + alpha_assoc * cv;
    if (event[i] == 1)
      target += weibull_lpdf(surv_time[i] | nu_weibull,
                             exp(-eta / nu_weibull));
    else
      target += weibull_lccdf(surv_time[i] | nu_weibull,
                              exp(-eta / nu_weibull));
  }

  // priors
  beta ~ normal(0, 10);
  sigma_b0 ~ exponential(1);
  sigma_b1 ~ exponential(1);
  sigma_y ~ exponential(1);
  alpha_trt ~ normal(0, 1);
  alpha_assoc ~ normal(0, 1);
  lambda_weibull ~ lognormal(0, 2);
  nu_weibull ~ gamma(2, 1);
}
