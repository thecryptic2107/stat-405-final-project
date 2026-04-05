//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] hba1c_tested;
  array[N] int readmitted_30;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real beta0;
  real beta1;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  beta0 ~ normal(0, 2);
  beta1 ~ normal(0, 1);
  readmitted_30 ~ bernoulli_logit(beta0 + beta1 * hba1c_tested);
}

generated quantities {
  vector[N] log_lik;
  array[N] int y_rep;

  for (i in 1:N) {
    real eta = beta0 + beta1 * hba1c_tested[i];
    log_lik[i] = bernoulli_logit_lpmf(readmitted_30[i] | eta);
    y_rep[i]   = bernoulli_logit_rng(eta);
  }
}

