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

data {
  int<lower=0> n;
  int<lower=0> p;
  vector[n] y;
  matrix[n,p] X;
}

parameters {
  vector[p] beta;
  real<lower=0> sig;
}

model {
  y ~ normal(X * beta, sig);
}
