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
  int<lower=0> N;
  int<lower=0> p;
  matrix[N,p] X;
  int<lower=0,upper=1> y[N];
}
parameters {
  vector[p] beta;
}
transformed parameters{
  vector[p] theta = inv_logit(beta);
}
model {
  y ~ bernoulli_logit(X * beta);
}

