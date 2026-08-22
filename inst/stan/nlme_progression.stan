// Nonlinear progression model (exponential decay shape)
data {
  int<lower=1> N;
  int<lower=1> n_obs;
  vector[n_obs] y;
  vector[n_obs] time;
  array[n_obs] int<lower=1, upper=N> id;
}
parameters {
  real b0;            // baseline
  real binf;          // asymptote
  real<lower=0> lambda;
  vector[N] z_b0;
  vector[N] z_slope;
  real<lower=0> sigma_b0;
  real<lower=0> sigma_slope;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[n_obs] yhat;
  for (j in 1:n_obs) {
    yhat[j] = (b0 + sigma_b0 * z_b0[id[j]]) +
              ((binf - b0) + sigma_slope * z_slope[id[j]]) *
              (1 - exp(-lambda * time[j]));
  }
}
model {
  y ~ normal(yhat, sigma_y);
  z_b0 ~ std_normal();
  z_slope ~ std_normal();
  b0 ~ normal(0, 20);
  binf ~ normal(0, 20);
  lambda ~ lognormal(0, 1);
  sigma_b0 ~ exponential(1);
  sigma_slope ~ exponential(0.5);
  sigma_y ~ exponential(1);
}
