# Fit the logistic regression described in P & C 
library(BayesLogit)  # For rpg function
library(mvtnorm)

#################
# Generate Data #
#################
# set.seed(071218)
n<-6
X<-matrix(c(0,1,1,0,0,
            1,1,0,1,0,
            1,1,0,0,1,
            1,0,1,1,0,
            1,0,1,0,1,
            0,0,0,1,1),
          ncol = 5,
          nrow = 6,
          byrow = TRUE)
p<-ncol(X)
y<-c(1,0,1,1,0,1)
fit<-glm(y~X,family=binomial)

# Priors
beta0<-rep(0,p)   # Prior mean of beta
T0<-diag(.01,p)   # Prior precision of beta

# Inits
beta<-rep(0,p)

#################
# Store Samples #
#################
nsim<-1000                # Number of MCMC Iterations
thin<-1				            # Thinning interval
burn<-nsim/2	            # Burnin
lastit<-(nsim-burn)/thin	# Last stored value
Beta<-matrix(0,lastit,p)

######### 
# Gibbs #
#########
tmp<-proc.time()          # Store current time

for (i in 1:nsim){
    eta<-X%*%beta
    w<-rpg(n,1,eta)
    z<-(y-1/2)/w                        # Or define z=y-1/2 and omit w in posterior mean m below
    v<-solve(crossprod(X*sqrt(w))+T0)   # Or solve(X%*%W%*%X), where W=diag(w) -- but this is slower
    m<-v%*%(T0%*%beta0+t(w*X)%*%z)      # Can omit w here if you define z=y-1/2
    beta<-c(rmvnorm(1,m,v))
    
    #################
    # Store Results #
    #################
    if (i> burn & i%%thin==0) {
        j<-(i-burn)/thin
        Beta[j,]<-beta 
    } 
    
    if (i%%100==0) print(i)
}

proc.time()-tmp             # MCMC run time

# Results
mbeta<-colMeans(Beta)
sbeta<-apply(Beta,2,sd)

# Compare MLEs and 
summary(fit)
cat("mbeta","\n",mbeta,"\n","sbeta","\n",sbeta)

# Trace plots
par(mfrow=c(2,1))
plot(1:lastit,Beta[,1],type="l",col="lightgreen")
abline(h=mbeta[1],col="blue4")

plot(1:lastit,Beta[,2],type="l",col="lightgreen")
abline(h=mbeta[2],col="blue4")

