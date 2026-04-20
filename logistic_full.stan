data {
  int<lower=0> N;
  vector[N] hba1c_tested;
  vector[N] age_middle;
  vector[N] age_old;
  vector[N] time_in_hospital;
  array[N] int<lower=0,upper=1> readmitted_30;
}

parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta3;
  real beta4;
}

model {
  // Priors — identical to Model 2 prior_sd <- c(2, 1, 1, 1, 1)
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  beta4 ~ normal(0, 1);

  // Likelihood via vectorised Bernoulli-logit
  readmitted_30 ~ bernoulli_logit(
    beta0 +
    beta1 * hba1c_tested +
    beta2 * age_middle +
    beta3 * age_old +
    beta4 * time_in_hospital
  );
}

generated quantities {
  // log_lik: used for LOO-CV / ELPD diagnostics
  vector[N] log_lik;
  // y_rep: posterior predictive draws for PPCs
  array[N] int y_rep;

  for (i in 1:N) {
    real eta = beta0
             + beta1 * hba1c_tested[i]
             + beta2 * age_middle[i]
             + beta3 * age_old[i]
             + beta4 * time_in_hospital[i];
    log_lik[i] = bernoulli_logit_lpmf(readmitted_30[i] | eta);
    y_rep[i]   = bernoulli_logit_rng(eta);
  }
}
