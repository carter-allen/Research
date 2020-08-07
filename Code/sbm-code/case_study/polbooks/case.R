library(sbmlogit)
library(igraph)
library(mcclust)
dyn.load("sbm-code/helper/newman.so")
dyn.load("sbm-code/helper/newman_ndc.so")
dyn.load("sbm-code/helper/nmi.so")

mp = function(vec, K){
  v = rep(1:K)
  l = length(vec)
  
  for (i in 1:K){
    v[i] = sum(vec==i)/l
  }
  
  return(v)
}

net <- 'polbooks'

# [ Parameters ]
g <- read.graph(paste0(net, '.gml'), format='gml')
nv <- as.integer(vcount(g)) 
ref <- scan(paste0('ref', net, '.txt'))

ng <- 10 # #(graphs)
tau2 <- 100
tol <- 1e-3
nsamples <- 1000
n0 <- 1

K <- as.integer(max(ref))
alpha <- rep(1/K, K)
runtime <- matrix(NA, nrow = ng, ncol = 3)
centroid.label <- matrix(NA, nrow = ng, ncol = nv)
binder.label <- matrix(NA, nrow = ng, ncol = nv)
kn.label <- matrix(NA, nrow = ng, ncol = nv)
centroid.sigma <- list()
centroid.gamma <- list()
centroid.eta <- list()

maprun <- 10
dir.create(paste0('results', maprun), showWarnings = F)

# sbm ============================
for(ir in 1:ng){
  message('=== ', ir, ' ===')
  lhood <- -Inf
  t.start <- proc.time()
  for(isig in 1:maprun){
    temp <- sbmlogit.map(g, alpha = K, infostep = 20)
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
  
  centroid.sigma[[ir]] <- res$sample
  centroid.gamma[[ir]] <- res$gamma
  centroid.eta[[ir]] <- res$eta
  
  # Newman w.o. degree correction ======================
  message('=== KN ===') 
  t.start <- proc.time() 
  snewman <- sbmlogit.remap(.Call('newman', g, K))
  t.end <- proc.time()
  kn.label[ir,] <- as.integer(snewman)
  runtime[ir, 3] <- as.numeric((t.end - t.start)[3])
  
  write.table(t(runtime[ir,]), file = paste0("results", maprun, "/runtime.txt"), 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  
  write.table(t(centroid.label[ir,]), paste0("results", maprun, "/centroid.txt"), 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  write.table(t(binder.label[ir,]), paste0("results", maprun, "/binder.txt") , 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T)
  write.table(t(kn.label[ir,]), paste0("results", maprun, "/kn.txt") , 
              row.names = FALSE,  col.names = FALSE, sep=" ", append = T) 
}
save(centroid.label, binder.label, kn.label, centroid.sigma, centroid.gamma, centroid.eta,
     file = paste0(net,maprun,".RData"))



ER <- function(x, y){
  return(sum(x != y) / length(x))
}

fiveNumber <- function(x){
  res <- matrix(NA, ncol(x), 5)
  for(i in 1:ncol(x)){
    res[i, ] <- quantile(x[, i], probs = c(.05, .25, .5, .75, .95))
  }
  return(res)
}

est <- list(centroid.label, binder.label, kn.label)
er.group.centroid <- matrix(NA, ng, K)
er.group.binder <- matrix(NA, ng, K)
er.group.kn <- matrix(NA, ng, K)
er <- matrix(NA, ng, 3)
nmi.group.centroid <- matrix(NA, ng, K)
nmi.group.binder <- matrix(NA, ng, K)
nmi.group.kn <- matrix(NA, ng, K)
nmi <- matrix(NA, ng, 3)
for(j in 1:ng){
  er[j, 1] <- ER(centroid.label[j, ], ref)
  er[j, 2] <- ER(binder.label[j, ], ref)
  er[j, 3] <- ER(kn.label[j, ], ref)
  nmi[j, 1] <- .Call("NMI0", as.integer(ref), as.integer(centroid.label[j, ]), nv, K)
  nmi[j, 2] <- .Call("NMI0", as.integer(ref), as.integer(binder.label[j, ]), nv, K)
  nmi[j, 3] <- .Call("NMI0", as.integer(ref), as.integer(kn.label[j, ]), nv, K)
}

for(i in 1:K){
  index <- which(ref == i)
  for(j in 1:ng){
    er.group.centroid[j, i] <- ER(centroid.label[j, index], ref[index])
    er.group.binder[j, i] <- ER(binder.label[j, index], ref[index])
    er.group.kn[j, i] <- ER(kn.label[j, index], ref[index])
    nmi.group.centroid[j, i] <- .Call("NMI0", as.integer(ref[index]), as.integer(centroid.label[j, index]), length(index), 1L)
    nmi.group.binder[j, i] <- .Call("NMI0", as.integer(ref[index]), as.integer(binder.label[j, index]), length(index), 1L)
    nmi.group.kn[j, i] <- .Call("NMI0", as.integer(ref[index]), as.integer(kn.label[j, index]), length(index), 1L)
  }
}


er
er.group.centroid
er.group.binder
er.group.kn

fiveNumber(er)
fiveNumber(er.group.centroid)
fiveNumber(er.group.binder)
fiveNumber(er.group.kn)
fiveNumber(nmi)
fiveNumber(nmi.group.centroid)
fiveNumber(nmi.group.binder)
fiveNumber(nmi.group.kn)
# save(er, er.group.binder, er.group.centroid, er.group.kn, file = 'errorRate.RData')


ests <- list(centroid.label, 
             binder.label, 
             kn.label, 
             centroid.woc.label, 
             centroid.ndc.label, 
             kn.ndc.label)

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

# save.image(file = '1.RData')


runtime <- read.table('results10/runtime.txt', header = F)
save(runtime, file = 'runtime.RData')


fiveNumber(runtime)
