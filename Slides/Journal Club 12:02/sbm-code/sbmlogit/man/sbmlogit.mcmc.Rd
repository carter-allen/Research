\name{sbmlogit.mcmc}
\alias{sbmlogit.mcmc}
\title{Gibbs sampler for Degree-Corrected Stochastic Blockmodels}
\description{
  Obtains posterior samples using a Gibbs sampler for degree-corrected
  stochastic blockmodels.
}
\usage{
sbmlogit.mcmc(graph, alpha = 2, tau2 = 100, size = 0, nsamples = 1,
infostep = 1, betastart = NULL, sigmastart = NULL, gammaconstraint = T)
}
\arguments{
  \item{graph}{
  an \code{igraph} graph object. If \code{size} is positive then each edge is
  assumed to have a binomial distribution where the number of trials is given
  in edge attribute \code{weight} in \code{graph}.
}
  \item{alpha}{
  a numeric vector of length \code{K} that specifies the hyper-parameters for
  the Dirichlet prior distribution on the probability of each node belonging
  to a community. For convenience, if \code{alpha} is a number then \code{K}
  is taken as that number and \code{alpha} is set to \code{rep(1, K)}.
}
  \item{tau2}{
  the prior variance on \code{beta}.
}
  \item{size}{
  the number of observations, if positive.
}
  \item{nsamples}{
  the number of samples to be generated, that is, the number of MCMC
  iterations.
}
  \item{infostep}{
  a background message is written every \code{infostep} steps to track
  convergence.
}
  \item{betastart}{
  initial value for \code{beta}. If not specified, \code{beta} will be
  randomly set according to its prior.
}
  \item{sigmastart}{
  initial value for \code{sigma}. If not specified, \code{sigma} will be
  randomly set according to its prior.
}
  \item{gammaconstraint}{
  logical specifying if \code{gamma} should be constrained for community
  detection or not. This is helpful when comparing to other models where the
  constraint is not set.
}
}
%\details{}

\value{
\code{sbmlogit.mcmc} returns a list with the following fields:

  \item{ngroups}{the number of communities, as specified by \code{alpha}.}
  \item{graph}{the input \code{graph}.}
  \item{pgroup}{the probabilities of belonging to each community.}
  \item{size}{the input \code{size}.}
  \item{gamma}{the MCMC samples for \code{gamma}.}
  \item{eta}{the MCMC samples for \code{eta}.}
  \item{sigma}{the MCMC samples for \code{sigma}.}
  \item{niters}{the number of iterations.}
  \item{lhood}{the values of the joint posterior at each iteration.}
  \item{mode}{the string \code{"mcmc"}. This field is used by
  \code{sbmlogit.plot}.}
}
\references{
For details on the specification of the degree-corrected stochastic
blockmodel, see "Bayesian Degree-Corrected Stochastic Blockmodels for
Community Detection".
}
%\author{}
%\note{}

\seealso{
\code{sbmlogit.map} for MAP estimation, and \code{sbmlogit.plot} to plot a
summary of the samples.
}

\examples{
library(igraphdata)

data(karate)
cm <- sbmlogit.map(karate, alpha=2, size=8)
cb <- sbmlogit.mcmc(karate, alpha=2, size=8, nsamples=1000)
\dontrun{sbmlogit.plot(cb, sbmlogit.remap(V(g))$Faction)}

data(UKfaculty)
g <- as.undirected(UKfaculty, mode="collapse")
cm <- sbmlogit.map(g, alpha=3)
cb <- sbmlogit.mcmc(g, alpha=3, nsamples=1000, sigmastart=cm$sigma)
\dontrun{sbmlogit.plot(cb, sbmlogit.remap(V(g))$Faction)}
}

\keyword{SBM}

