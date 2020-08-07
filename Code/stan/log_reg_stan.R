library(rstan)
# rstan code for fitting logistic regression
# relies on log_reg.stan

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/Documents/School/Research/Code/stan")

# simulate logistic regression data
n = 500 # number of observations
p = 3 # number of predictors
beta = c(1, -2, 1.5) # regression parameters
X = matrix(rnorm(n*p), nrow = n, ncol = p) # predictor matrix
X[,1] = 1 # force an intercept
mu = X %*% beta # logit(eta)
eta = exp(mu)/(1 + exp(mu)) # expit(mu)
y = rbinom(n,1,eta) # random obs

# define stan data list
stan_dat <- list(n = n,
                 X = X,
                 y = y,
                 p = p)

# fit regression model
fit <- stan(file = "log_reg.stan",
            data = stan_dat)