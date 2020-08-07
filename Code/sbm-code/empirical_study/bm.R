library(sbmlogit)
library(igraph)
library(mcclust)
dyn.load("../helper/newman.so")

args <- (commandArgs(TRUE))

# [ Parameters ]
# nv <- as.integer(100) # #(nodes): 100, 500
# adeg <- 10 # average degree: 10, 15, 25
# gamma <- 2 # exponent of degree distribution: 2, 3
# beta <- 1 # exponent of community size distribution: 1, 2
nv <- as.integer(args[1]) # #(nodes): 100, 500
adeg <- args[2] # average degree: 10, 15, 25
gamma <- args[3] # exponent of degree distribution: 2, 3
beta <- args[4] # exponent of community size distribution: 1, 2
mdeg <- as.integer(nv / 2) # max degree
mu <- (1:6) / 10 # mixing parameter: .1, .2, ..., .5 [*],  .6

ng <- 100 # #(graphs)
tau2 <- 100
tol <- 1e-3
nsamples <- 1000
n0 <- 1

mp = function(vec, K){
  v = rep(1:K)
  l = length(vec)
  
  for (i in 1:K){
    v[i] = sum(vec==i)/l
  }
  
  return(v)
}

# [ Benchmark ]
inv <- 1
runtime <- matrix(0, ng*length(mu), 5)

for(ib in 1:length(beta)){
  for (ig in 1:length(gamma)) {
    for (ia in 1:length(adeg)) {
      dirName <- paste0('nv', nv[inv], 'a', gamma[ig], 'b', beta[ib], 'deg', adeg[ia])
      dir.create(paste0(dirName, '/results'), showWarnings = F)
      load(paste0(dirName, '/reflabel.RData'))
      centroid.woc.label <- matrix(NA, ng*length(mu), nv[inv])
      binder.woc.label <- matrix(NA, ng*length(mu), nv[inv])
      centroid.ndc.label <- matrix(NA, ng*length(mu), nv[inv])
      binder.ndc.label <- matrix(NA, ng*length(mu), nv[inv])
      kn.ndc.label <- matrix(NA, ng*length(mu), nv[inv])
      
      for (im in 1:length(mu)) {
        for (ir in 1:ng) {
          cat(paste0("\n=== beta = ", beta[ib], ", gamma = ", gamma[ig], ", adeg = ", adeg[ia],
                     ", mu = ", mu[im], ", graph = ", ir, " ===\n"))
          
          ref <- ref.label[ng*(im-1) + ir,]
          K <- as.integer(max(ref))
          alpha <- rep(1/K, K)
          g <- read.graph(paste0(dirName, '/graph/mu', mu[im], 'graph', ir, '.gml'), format='gml')
          
          # without constraint ============================
          message('=== woc ===')  
          lhood <- -Inf
          t.start <- proc.time()
          for(isig in 1:20){
            temp <- sbmlogit.map(g, alpha = K, infostep = 5, gammaconstraint = F)
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
          centroid.woc.label[(ng*(im-1)+ir),] <- as.integer(sbmlogit.remap(scentroid))
          runtime[(ng*(im-1)+ir), 1] <- as.numeric((t.end - t.start)[3]) + as.numeric((t.end1 - t.start1)[3])
          
          # binder -------------
          # message('start binder ...')
          t.start1 <- proc.time()
          sbinder <- minbinder(comp.psm(as.matrix(Sigma)), "laugreen")$cl
          t.end1 <- proc.time()
          binder.woc.label[(ng*(im-1)+ir),] <- as.integer(sbmlogit.remap(sbinder))
          runtime[(ng*(im-1)+ir), 2] <- as.numeric((t.end - t.start)[3]) + as.numeric((t.end1 - t.start1)[3])
          
          # without degree correction ======================
          message('=== ndc ===') 
          lhood <- -Inf
          t.start <- proc.time()
          for(isig in 1:20){
            temp <- sbmlogit.map.ndc(g, alpha = K, infostep = 5)
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
          centroid.ndc.label[(ng*(im-1)+ir),] <- as.integer(sbmlogit.remap(scentroid))
          runtime[(ng*(im-1)+ir), 3] <- as.numeric((t.end - t.start)[3]) + as.numeric((t.end1 - t.start1)[3])
          
          # binder -------------
          # message('start binder ...')
          t.start1 <- proc.time()
          sbinder <- minbinder(comp.psm(as.matrix(Sigma)), "laugreen")$cl
          t.end1 <- proc.time()
          binder.ndc.label[(ng*(im-1)+ir),] <- as.integer(sbmlogit.remap(sbinder))
          runtime[(ng*(im-1)+ir), 4] <- as.numeric((t.end - t.start)[3]) + as.numeric((t.end1 - t.start1)[3])
          
          # Newman w.o. degree correction ======================
          message('=== KNndc ===') 
          t.start <- proc.time() 
          snewman <- sbmlogit.remap(.Call('newman', g, K))
          t.end <- proc.time()
          kn.ndc.label[(ng*(im-1)+ir),] <- as.integer(snewman)
          runtime[(ng*(im-1)+ir), 5] <- as.numeric((t.end - t.start)[3])
          
          write.table(t(runtime[(ng*(im-1)+ir),]), paste0(dirName,"/results/runtimeMcmcOther.txt") , 
                      row.names = F,  col.names = F, sep=" ", append = T)
          
          write.table(t(centroid.woc.label[(ng*(im-1)+ir),]), paste0(dirName,"/results/centroidWoc.txt") , 
                      row.names = F,  col.names = F, sep=" ", append = T)
          write.table(t(binder.woc.label[(ng*(im-1)+ir),]), paste0(dirName,"/results/binderWoc.txt") , 
                      row.names = F,  col.names = F, sep=" ", append = T)
          
          write.table(t(centroid.ndc.label[(ng*(im-1)+ir),]), paste0(dirName,"/results/centroidNdc.txt") , 
                      row.names = F,  col.names = F, sep=" ", append = T)
          write.table(t(binder.ndc.label[(ng*(im-1)+ir),]), paste0(dirName,"/results/binderNdc.txt") , 
                      row.names = F,  col.names = F, sep=" ", append = T)
          
          write.table(t(kn.ndc.label[(ng*(im-1)+ir),]), paste0(dirName,"/results/knNdc.txt") , 
                      row.names = F,  col.names = F, sep=" ", append = T)  
        }
      }
      save(runtime, centroid.woc.label, binder.woc.label, centroid.ndc.label, binder.ndc.label, kn.ndc.label, 
           file = paste0(dirName,"/results/mcmcOther.RData"))
    }
  }
}
