# run stan analysis example
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#how-to-use-rstan
setwd("~/Documents/School/Research/Code/stan")

library("rstan") # observe startup messages
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# define data as list
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

# fit model
fit <- stan(file = 'schools.stan', data = schools_dat)

# plot fit
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))

la <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- la$mu 

### return an array of three dimensions: iterations, chains, parameters 
a <- extract(fit, permuted = FALSE) 

### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)
