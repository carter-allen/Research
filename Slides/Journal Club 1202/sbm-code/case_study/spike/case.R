library(sbmlogit)
library(igraph)
library(mcclust)
dyn.load("../../helper/newman.so")
dyn.load("../../helper/newman_ndc.so")
dyn.load("../../helper/nmi.so")

mp = function(vec, K){
  v = rep(1:K)
  l = length(vec)
  
  for (i in 1:K){
    v[i] = sum(vec==i)/l
  }
  
  return(v)
}

net <- 'spike'

# [ Parameters ]
g <- read.graph(paste0(net, '.gml'), format='gml')
nv <- as.integer(vcount(g)) 
ref <- scan(paste0('ref', net, '.txt'))

ng <- 100 # #(graphs)
tau2 <- 100
tol <- 1e-3
nsamples <- 1000
n0 <- 1

K <- as.integer(max(ref))
alpha <- rep(1/K, K)
runtime <- matrix(NA, nrow = ng, ncol = 6)
centroid.label <- matrix(NA, nrow = ng, ncol = nv)
binder.label <- matrix(NA, nrow = ng, ncol = nv)
kn.label <- matrix(NA, nrow = ng, ncol = nv)
centroid.woc.label <- matrix(NA, nrow = ng, ncol = nv)
centroid.ndc.label <- matrix(NA, nrow = ng, ncol = nv)
kn.ndc.label <- matrix(NA, nrow = ng, ncol = nv)
centroid.gamma <- list()
centroid.eta <- list()
centroid.woc.gamma <- list()
centroid.woc.eta <- list()
centroid.ndc.gamma <- list()
centroid.ndc.eta <- list()

dir.create('results1', showWarnings = F)

# sbm ============================
for(ir in 1:ng){
  message('=== ', ir,' ===')
  lhood <- -Inf
  t.start <- proc.time()
  for(isig in 1:20){
    temp <- sbmlogit.map(g, alpha = K, infostep = 10)
    lhood.temp <- temp$lhood[length(temp$lhood)]
    if (lhood.temp > lhood){
      sigma0 <- temp$sigma
      lhood <- lhood.temp
    }
  }
  res <- sbmlogit.mcmc(graph = g, alpha = K, tau2 = tau2, 
                       nsamples = nsamples, infostep = 100,
                       sigmastart = sigma0)
  t.end <- proc.time()
  conv <- floor(.2*nsamples) + 1
  Sigma <- res$sample[conv : nsamples, ]
  
  # centroid -----------
  # message('start centroid ...')          
  t.start1 <- proc.time()
  scentroid <- apply(t(apply(Sigma, 2, mp, K)), 1, which.max)
  t.end1 <- proc.time()
  centroid.label[ir,] <- as.integer(sbmlogit.remap(scentroid))
  runtime[ir, 1] <- as.numeric((t.end - t.start)[3]) + as.numeric((t.end1 - t.start1)[3])
  
  # binder -------------
  # message('start binder ...')
  t.start1 <- proc.time()
  sbinder <- minbinder(comp.psm(as.matrix(Sigma)), "laugreen")$cl
  t.end1 <- proc.time()
  binder.label[ir,] <- as.integer(sbmlogit.remap(sbinder))
  runtime[ir, 2] <- as.numeric((t.end - t.start)[3]) + as.numeric((t.end1 - t.start1)[3])
  
  #   centroid.sigma[[ir]] <- res$sample
  centroid.gamma[[ir]] <- res$gamma
  centroid.eta[[ir]] <- res$eta
  
  # Newman w degree correction ======================
  message('=== KN ===') 
  t.start <- proc.time() 
  snewman <- sbmlogit.remap(.Call('newman', g, K))
  t.end <- proc.time()
  kn.label[ir,] <- as.integer(snewman)
  runtime[ir, 3] <- as.numeric((t.end - t.start)[3])
  
  
  # without constraint ============================
  message('=== woc ===')  
  lhood <- -Inf
  t.start <- proc.time()
  for(isig in 1:20){
    temp <- sbmlogit.map(g, alpha = K, infostep = 10, gammaconstraint = F)
    lhood.temp <- temp$lhood[length(temp$lhood)]
    if (lhood.temp > lhood){
      sigma0 <- temp$sigma
      lhood <- lhood.temp
    }
  }                       
  res <- sbmlogit.mcmc(graph = g, alpha = K, tau2 = tau2, 
                       nsamples = nsamples, infostep = 100,
                       sigmastart = sigma0, gammaconstraint = F)
  t.end <- proc.time()
  conv = floor(.2*nsamples) + 1
  Sigma = res$sample[conv : nsamples, ]
  
  # centroid -----------
  # message('start centroid ...')          
  t.start1 <- proc.time()
  scentroid <- apply(t(apply(Sigma, 2, mp, K)), 1, which.max)
  t.end1 <- proc.time()
  centroid.woc.label[ir,] <- as.integer(sbmlogit.remap(scentroid))
  runtime[ir, 4] <- as.numeric((t.end - t.start)[3]) + as.numeric((t.end1 - t.start1)[3])
  
  #   centroid.woc.sigma[[ir]] <- res$sample
  centroid.woc.gamma[[ir]] <- res$gamma
  centroid.woc.eta[[ir]] <- res$eta
  
  # without degree correction ======================
  message('=== ndc ===') 
  lhood <- -Inf
  t.start <- proc.time()
  for(isig in 1:20){
    temp <- sbmlogit.map.ndc(g, alpha = K, infostep = 10)
    lhood.temp <- temp$lhood[length(temp$lhood)]
    if (lhood.temp > lhood){
      sigma0 <- temp$sigma
      lhood <- lhood.temp
    }
  }
  res <- sbmlogit.mcmc.ndc(graph = g, alpha = K, tau2 = tau2, 
                           nsamples = nsamples, infostep = 100,
                           sigmastart = sigma0)
  t.end <- proc.time()
  conv = floor(.2*nsamples) + 1
  Sigma = res$sample[conv : nsamples, ]
  
  # centroid -----------
  # message('start centroid ...')
  t.start1 <- proc.time()
  scentroid <- apply(t(apply(Sigma, 2, mp, K)), 1, which.max)
  t.end1 <- proc.time()
  centroid.ndc.label[ir,] <- as.integer(sbmlogit.remap(scentroid))
  runtime[ir, 5] <- as.numeric((t.end - t.start)[3]) + as.numeric((t.end1 - t.start1)[3])
  
  #   centroid.sigma[[ir]] <- res$sample
  centroid.ndc.gamma[[ir]] <- res$gamma
  centroid.ndc.eta[[ir]] <- res$eta
  
  # Newman w degree correction ======================
  message('=== KN ===') 
  t.start <- proc.time() 
  snewman <- sbmlogit.remap(.Call('newman1', g, K))
  t.end <- proc.time()
  kn.ndc.label[ir,] <- as.integer(snewman)
  runtime[ir, 6] <- as.numeric((t.end - t.start)[3])
  
  
  write.table(t(runtime[ir,]), file = paste0("results1/runtime.txt"), 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  
  write.table(t(centroid.label[ir,]), paste0("results1/centroid.txt"), 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  write.table(t(binder.label[ir,]), paste0("results1/binder.txt") , 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  write.table(t(kn.label[ir,]), paste0("results1/kn.txt") , 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  
  write.table(t(centroid.woc.label[ir,]), paste0("results1/centroidWoc.txt"), 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  write.table(t(centroid.ndc.label[ir,]), paste0("results1/centroidNdc.txt"), 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  write.table(t(kn.ndc.label[ir,]), paste0("results1/knNdc.txt"), 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  
}


save(centroid.label, binder.label, kn.label, 
     centroid.woc.label, centroid.ndc.label, kn.ndc.label, 
     centroid.gamma, centroid.eta,
     centroid.woc.gamma, centroid.woc.eta,
     centroid.ndc.gamma, centroid.ndc.eta,
     file = paste0(net,"1.RData"))

fiveNumber <- function(x){
  res <- matrix(NA, ncol(x), 5)
  for(i in 1:ncol(x)){
    res[i, ] <- quantile(x[, i], probs = c(.05, .25, .5, .75, .95))
  }
  return(res)
}

gamma <- matrix(NA, ng, choose(K,2))
gamma.woc <- matrix(NA, ng, choose(K,2))
gamma.ndc <- matrix(NA, ng, choose(K,2))
for(j in 1:ng){
  gamma[j, ] <- mean(centroid.gamma[[j]][201:1000])
  gamma.woc[j, ] <- mean(centroid.woc.gamma[[j]][201:1000])
  gamma.ndc[j, ] <- mean(centroid.ndc.gamma[[j]][201:1000])
}
fiveNumber(gamma)
fiveNumber(gamma.woc)
fiveNumber(gamma.ndc)

eta <- matrix(NA, ng, nv)
eta.woc <- matrix(NA, ng, nv)
eta.ndc <- matrix(NA, ng, 1)
for(j in 1:ng){
  eta[j, ] <- colMeans(centroid.eta[[j]][201:1000,])
  eta.woc[j, ] <- colMeans(centroid.woc.eta[[j]][201:1000,])
  eta.ndc[j, ] <- mean(centroid.ndc.eta[[j]][201:1000])
}
fiveNumber(eta)
fiveNumber(eta.woc)
fiveNumber(eta.ndc)


ests <- list(centroid.label, binder.label, kn.label, 
             centroid.woc.label, centroid.ndc.label, kn.ndc.label)

er <- matrix(NA, ng, length(ests))
nmi <- matrix(NA, ng, length(ests))

for(i in 1:length(ests)){
  for(j in 1:ng){
    er[j, i] <- mean(ests[[i]][j, ] != ref)
    nmi[j,i] <- .Call("NMI0", as.integer(ref), as.integer(sbmlogit.remap(ests[[i]][j, ])), nv, K)
  }
}

fiveNumber(er)
fiveNumber(nmi)
