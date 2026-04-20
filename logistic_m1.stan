data {
  int<lower=0> N;
  vector[N] hba1c_tested;
  array[N] int<lower=0,upper=1> readmitted_30;
}

parameters {
  real beta0;
  real beta1;
}

model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  readmitted_30 ~ bernoulli_logit(beta0 + beta1 * hba1c_tested);
}
