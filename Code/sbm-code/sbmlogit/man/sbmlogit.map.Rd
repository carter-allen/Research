\name{sbmlogit.map}
\alias{sbmlogit.map}
\title{MAP Estimation for Degree-Corrected Stochastic Blockmodels}
\description{
  Obtains a MAP estimate for degree-corrected stochastic blockmodels.
}
\usage{
sbmlogit.map(graph, alpha = 2, tau2 = 100, size = 0, tolerance = 1e-06,
maxiter = 25, infostep = 1, betastart = NULL, sigmastart = NULL,
gammaconstraint = T)
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
  \item{tolerance}{
  positive tolerance for assessing convergence. The method converges when
  \code{|logprob - logprob_{old}| > tolerance}, where \code{logprob} is the
  log of the joint posterior distribution.
}
  \item{maxiter}{
  maximum number of iterations.
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
\code{sbmlogit.map} returns a list with the following fields:

  \item{ngroups}{the number of communities, as specified by \code{alpha}.}
  \item{graph}{the input \code{graph}.}
  \item{pgroup}{the probabilities of belonging to each community.}
  \item{size}{the input \code{size}.}
  \item{gamma}{the mode estimate for \code{gamma}.}
  \item{eta}{the mode estimate for \code{eta}.}
  \item{sigma}{the mode estimate for \code{sigma}.}
  \item{sd}{a local estimate for the standard deviations of \code{gamma} and
  \code{eta} from the iterative weighted least squares.}
  \item{niters}{the number of iterations until convergence.}
  \item{delta}{the last value for \code{|logprob - logprob_{old}|}, which can
  be used to say if the method converged or if \code{maxiter} was reached.}
  \item{lhood}{the values of the joint posterior at each iteration.}
  \item{mode}{the string \code{"map"}. This field is used by
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
\code{sbmlogit.mcmc} for posterior sampling using MCMC, and
\code{sbmlogit.plot} to plot the fit.
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

