# Community detection through GLM
library(igraph)

# [ Parameters ]
EPSILON <- 1e-6 # regularization

SWEEP <- function (n, m, A, ind) # A is double, n-by-m
  .Call(S_sweep, as.integer(n), as.integer(m), A, as.integer(ind - 1))

rtruncnorm <- function (mu, sigma)
  .Call(S_rtruncnorm, as.double(mu), as.double(sigma))


rdirichlet <- function (alpha) {
  s <- rgamma(length(alpha), alpha, 1)
  s / sum(s)
}

sample.sigma <- function (n, K) {
  repeat {
    sigma <- sapply(rep(K, n), sample.int, 1)
    cpi <- tapply(sigma, factor(sigma, 1:K), length)
    if (all(!is.na(cpi) & cpi > 1)) break
  }
  as.integer(sigma - 1)
}


# [ update based on *active set method* for gamma < 0 ]
update.beta <- function (d, gamma.ind, active, beta, m, S) {
  lambda <- m[active] # Lagrange multipliers
  if (any(lambda < 0)) {
    j <- which(lambda < 0)
    active <- active[-j] # remove constraint
  }

  # compute step:
  non.active <- which(is.na(match(1:d, active)))
  C <- chol(matrix(S, nrow = d, ncol = d))
  V <- matrix(SWEEP(d, d, as.double(chol2inv(C)), active), nrow = d)
  if (any(is.na(V)))
    V[is.na(V)] <- 0 # FIXME: workaround
  p <- V[non.active, non.active] %*% m[non.active] # step
  if (any(is.na(p)))
    stop("invalid step")
  V <- diag(V); V[active] <- 0

  alpha.min <- 1
  cand <- which(non.active %in% gamma.ind & p > 0)
  if (length(cand) > 0) {
    offset <- -EPSILON
    alpha <- (offset - beta[non.active[cand]]) / p[cand]
    alpha.min <- min(1, alpha)
    if (alpha.min < 1) { # any blocking constraints?
      j <- which(alpha == alpha.min)
      active <- c(active, non.active[cand[j]]) # activate constraint
    }
  }
  if (is.na(alpha.min))
    stop("invalid min alpha")
  beta[non.active] <- beta[non.active] + alpha.min * p

  return(list(beta=beta, active=active, sd=sqrt(V)))
}


sbmlogit.remap <- function (sigma) {
  k <- length(sigma)
  visited <- rep(F, k)
  revind <- rep(-1, k)
  c <- 1
  for (i in 1:k) {
    l <- sigma[i]
    if (!visited[l]) { # first visit?
      revind[l] <- c # order of appearance
      c <- c + 1
      visited[l] <- T
    }
  }
  for (l in 1:k) {
    if (revind[l] == -1) { # no appearance?
      revind[l] <- c
      c <- c + 1
    }
  }
  for (i in 1:length(sigma))
    sigma[i] <- revind[sigma[i]]
  return(sigma)
}


# [ Auxiliary ]
colprob <- function (col, p) {
  r <- c(col2rgb(col) / 255)
  r <- (1 - p) * c(1, 1, 1) + p * r
  r <- ifelse(r < 0, 0, r) # truncate errors
  rgb(r[1], r[2], r[3])
}

darken <- function (col, p = .5) {
  r <- p * c(col2rgb(col) / 255)
  rgb(r[1], r[2], r[3])
}

circle <- function (center, radius, col, border, lwd=2, n=100) {
  t <- seq(0, 2 * pi, length.out = n)
  x <- center[1] + radius * cos(t)
  y <- center[2] + radius * sin(t)
  polygon(x, y, col=col, border=border, lwd=lwd)
}

logit <- function (x) log(x / (1 - x))
ilogit <- function (x) 1 / (1 + exp(-x))

# `k` is number of groups
# `p` is data frame with marginal posterior for each group
rfactor <- 300
plotgraph <- function (g, p, ref, xg, lwd=1) {
  # graph
  if (missing(xg))
    xg <- layout.fruchterman.reingold(g) # coordinates
  minr <- max(apply(xg, 2, max) - apply(xg, 2, min)) / rfactor # min radius
  maxr <- 10 * minr # max radius
  dg <- degree(g)
  sg <- minr + (maxr - minr) / (max(dg) - min(dg)) * dg # scale by degree
  # vertex colors
  label <- apply(p, 1, which.max)
  k <- dim(p)[2]
  rp <- (apply(p, 1, max) - 1 / k) / (1 - 1 / k) # relative proportions
  # plot
  plot(xg, type="n", axes=F, xlab="", ylab="", frame=F, asp=1)
  eg <- get.edgelist(g)
  for (i in 1:dim(eg)[1]) {
    e1 <- eg[i, 1]; e2 <- eg[i, 2]
    lines(c(xg[e1, 1], xg[e2, 1]), c(xg[e1, 2], xg[e2, 2]), col="gray")
  }
  label <- label + 1 # color code
  if (missing(ref)) {
    for (i in 1:length(V(g)))
      circle(xg[i,], sg[i], colprob(label[i], rp[i]), "black", lwd=1)
  }
  else {
    ref <- ref + 1 # color code
    for (i in 1:length(V(g)))
      circle(xg[i,], sg[i], colprob(label[i], rp[i]), darken(ref[i]), lwd)
  }
  invisible(xg)
}

sbmlogit.plot <- function (cm, ref=NULL) {
  ngroups <- cm$ngroups
  g <- cm$graph
  # compute p
  p <- matrix(nrow = vcount(g), ncol = ngroups)
  if (cm$mode == "map") {
    for (k in 1:ngroups)
      p[,k] <- ifelse(cm$sigma == k, 1, 0)
  }
  else {
    for (k in 1:ngroups)
      p[,k] <- apply(cm$sample, 2, function(x) sum(x==k) / cm$niters)
  }
  # compute eta and label (estimate)
  if (cm$mode == "mcmc") {
    eta <- apply(cm$eta, 2, mean)
    label <- apply(p, 1, which.max)
  }
  else {
    eta <- cm$eta
    label <- cm$sigma
  }
  
  # report
  if (!is.null(ref))
    print(c("Recall", sum(label==ref) / length(label)))

  if (cm$mode == "mcmc") {
    op <- par(mfrow=c(1, 3))
    # [1]
    plot(apply(p, 1, max), col=label + 1, pch=19, ylim=c(1 / ngroups, 1),
         xlab='Vertex', ylab='Posterior')
    abline(h=.5, lty=2)
  }
  else
    op <- par(mfrow=c(1, 2))
  # [2]
  xg <- plotgraph(g, p, ref,, 2)
  # [3]
  if (cm$size == 0) {
    plot(logit(degree(g) / (vcount(g) - 1)), eta, col=label + 1, pch=19,
         xlab="logit(degree)", ylab="eta")
  } else {
    value <- E(g)$value # FIXME
    s <- graph.strength(g, weights=value)
    plot(logit(s / (sum(s) - s)), eta, col=label + 1,
         pch=19, xlab="logit(weighted degree)", ylab="eta")
  }
  par(op)
  invisible(list(x=xg, p=p))
}



# [ Main ]

sbmlogit.nc.map <- function (graph, tau2 = 100, size = 0,
                          tolerance = 1e-6, maxiter = 25, infostep = 1,
                          betastart = NULL) {
  # [ init ]
  n <- vcount(graph)
  K <- as.integer(1) # no community
  tau2 <- as.double(tau2)
  if (is.null(betastart))
    beta <- numeric(n)
  else
    beta <- betastart
  S <- as.double(matrix(nrow = n, ncol = n))
  size <- as.double(size)

  # [ fit ]
  is <- 0
  sample <- as.integer(F) # MAP
  lhood <- numeric(maxiter)
  lcur <- -Inf
  repeat {
    is <- is + 1
    # [ fit beta ]
    last.beta <- beta
    m <- .Call(S_comm_beta_map_bin, graph, NULL, K, tau2,
               as.double(beta), size, S)
    C <- chol(matrix(S, nrow = n, ncol = n))
    m <- backsolve(C, backsolve(C, m, transpose=T)) # m <- solve(S) * m
    beta <- beta + m
    sd <- diag(C)
    lhood[is] <- - .5 * (sum(beta ^ 2) / tau2)

    # [ report ]
    if ((infostep > 0) && (is %% infostep == 0))
      message("[", is, "] lhood = ", lhood[is])

    delta <- abs(lcur - lhood[is])
    if (is > maxiter || delta < tolerance) break
    lcur <- lhood[is]
  }
  return(list(ngroups=K, graph=graph, size=size, eta=last.beta,
              sd=sd, niters=is, delta=delta, lhood=lhood[1:is], mode="map"))
}


sbmlogit.map <- function (graph, alpha = 2, tau2 = 100, size = 0,
                          tolerance = 1e-6, maxiter = 25, infostep = 1,
                          betastart = NULL, sigmastart = NULL,
                          gammaconstraint = T) {
  # [ init ]
  n <- vcount(graph)
  if (length(alpha) == 1) # alpha is K?
    alpha <- rep(1, alpha)
  logpi <- as.double(log(alpha / sum(alpha))) # log(prior mean)
  K <- as.integer(length(logpi))
  Ck <- choose(K, 2)
  tau2 <- as.double(tau2)
  if (is.null(betastart))
    # avoid binding constraints:
    beta <- rep(-EPSILON, n + Ck)
  else
    beta <- betastart
  if (is.null(sigmastart))
    sigma <- sample.sigma(n, K)
  else
    sigma <- as.integer(sigmastart - 1)
  S <- as.double(matrix(nrow = n + Ck, ncol = n + Ck))
  size <- as.double(size)
  gamma.ind <- 1:Ck
  eta.ind <- Ck + 1:n
  active <- which(beta[gamma.ind] >= 0)

  # [ fit ]
  is <- 0
  sample <- as.integer(F) # MAP
  lhood <- numeric(maxiter)
  lcur <- -Inf
  repeat {
    is <- is + 1
    # [ fit beta ]
    last.beta <- beta
    m <- .Call(S_comm_beta_map_bin, graph, sigma, K, tau2,
               as.double(beta), size, S)
    if (gammaconstraint) {
      u <- update.beta(n + Ck, gamma.ind, active, beta, m, S)
      beta <- u$beta; active <- u$active; sd <- u$sd
    }
    else {
      C <- chol(matrix(S, nrow = n + Ck, ncol = n + Ck))
      m <- backsolve(C, backsolve(C, m, transpose=T)) # m <- solve(S) * m
      beta <- beta + m
      sd <- diag(C)
    }

    # [ fit sigma ]
    last.sigma <- sigma
    lhood[is] <- .Call(S_comm_sigma_bin, graph, sigma, logpi,
                       as.double(beta), size, sample)
    lhood[is] <- lhood[is] - .5 * (sum(beta ^ 2) / tau2)

    # [ fit pi ]
    last.logpi <- logpi
    cpi <- tapply(sigma, factor(sigma, 0:(K-1)), length)
    cpi[is.na(cpi)] <- 1
    cpi <- cpi + alpha - 1  # counts
    logpi <- as.double(log(cpi / sum(cpi)))
    lhood[is] <- lhood[is] + sum((alpha - 1) * logpi)

    # [ report ]
    if ((infostep > 0) && (is %% infostep == 0))
      message("[", is, "] lhood = ", lhood[is])

    delta <- abs(lcur - lhood[is])
    if (is > maxiter || delta < tolerance) break
    lcur <- lhood[is]
  }
  return(list(ngroups=K, graph=graph, pgroup=exp(logpi), size=size,
              gamma=last.beta[gamma.ind], eta=last.beta[eta.ind],
              sigma=last.sigma + 1, sd=sd, niters=is, delta=delta,
              lhood=lhood[1:is], mode="map"))
}


sbmlogit.mcmc <- function (graph, alpha = 2, tau2 = 100, size = 0,
                           nsamples = 1, infostep = 1,
                           betastart = NULL, sigmastart = NULL,
                           gammaconstraint = T) {
  # [ init ]
  n <- vcount(graph)
  if (length(alpha) == 1) # alpha is K?
    alpha <- rep(1, alpha)
  logpi <- as.double(log(alpha / sum(alpha))) # log(prior mean)
  K <- as.integer(length(logpi))
  Ck <- choose(K, 2)
  tau2 <- as.double(tau2)
  if (is.null(betastart))
    beta <- rep(0, n + Ck)
  else
    beta <- betastart
  if (is.null(sigmastart))
    sigma <- as.integer(sapply(rep(K, n), sample.int, 1) - 1)
  else
    sigma <- as.integer(sigmastart - 1)

  S <- as.double(matrix(nrow = n + Ck, ncol = n + Ck))
  size <- as.double(size)
  gamma.ind <- 1:Ck
  eta.ind <- Ck + 1:n

  # [ sample ]
  sample <- as.integer(T)
  ss <- matrix(nrow = nsamples, ncol = n)
  bs <- matrix(nrow = nsamples, ncol = n + Ck)
  ps <- matrix(nrow = nsamples, ncol = K)
  lhood <- numeric(nsamples)
  for (is in 1:nsamples) {
    # [ sample beta ]
    m <- .Call(S_comm_beta_mcmc_bin, graph, sigma, K, tau2,
               as.double(beta), size, S)
    C <- chol(matrix(S, nrow = n + Ck, ncol = n + Ck))
    beta <- backsolve(C, backsolve(C, m, transpose=T))

    # sample eta marginally
    delta <- backsolve(C[eta.ind, eta.ind], rnorm(length(eta.ind)))
    beta[eta.ind] <- beta[eta.ind] + delta
    # sample gamma | eta
    C.invdiag <- as.double(1 / diag(C)[gamma.ind])
    beta.mu <- beta[gamma.ind] + (C[gamma.ind, eta.ind] %*% delta) * C.invdiag
    if (gammaconstraint)
      beta[gamma.ind] <- -.Call(S_rtruncnorm, as.double(-beta.mu), C.invdiag)
    else
      beta[gamma.ind] <- -rnorm(length(gamma.ind), -beta.mu, C.invdiag)
    bs[is,] <- beta

    # [ sample sigma ]
    lhood[is] <- .Call(S_comm_sigma_bin, graph, sigma, logpi,
                       as.double(beta), size, sample)
    lhood[is] <- lhood[is] - .5 * (sum(beta ^ 2) / tau2)
    ss[is,] <- sigma + 1

    # [ sample pi ]
    cpi <- tapply(sigma, factor(sigma, 0:(K-1)), length)
    cpi[is.na(cpi)] <- 1
    ppi <- rdirichlet(cpi + alpha)  # counts
    logpi <- as.double(log(ppi))
    lhood[is] <- lhood[is] + sum((alpha - 1) * logpi)
    ps[is,] <- ppi

    # [ report ]
    if (is %% infostep == 0)
      message("[", is, "] lhood = ", lhood[is])
  }
  return(list(ngroups=K, graph=graph, niters=is, size=size,
              gamma=bs[,gamma.ind], eta=bs[,eta.ind], pgroup=ps,
              sample=ss, lhood=lhood, mode="mcmc"))
}

