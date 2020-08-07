library(rstan)
# rstan code for fitting linear regression
# relies on lin_reg.stan

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/Documents/School/Research/Code/stan")

# simulate linear regression data
n = 1000 # number of observations
p = 3 # number of predictors (including intercept)
beta = c(1,-4,10) # true coefs
sig = 2 # sd of resids
X = matrix(rnorm(n*p), nrow = n, ncol = p) # predictor matrix
X[,1] = 1 # force an intercept
e = rnorm(n,0,sig) # residuals
y = X %*% beta + e # simulated responses
y = as.vector(y)

# define stan data list
stan_dat <- list(n = n,
    X = X,
    y = y,
    p = p)

# fit regression model
fit <- stan(file = "lin_reg.stan",
            data = stan_dat)
