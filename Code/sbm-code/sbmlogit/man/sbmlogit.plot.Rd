\name{sbmlogit.plot}
\alias{sbmlogit.plot}
\title{Plotting stochastic blockmodels and graphs}
\description{
  Given a stochastic blockmodel specification on a graph, plots the graph with
  community labels and a summary of the parameters.
}
\usage{sbmlogit.plot(cm, ref = NULL)}
\arguments{
  \item{cm}{
  a list as returned by \code{sbmlogit.map} or \code{sbmlogit.mcmc}.
}
  \item{ref}{
  a label configuration vector that specifies a reference. The border of
  vertices in the graph plot are colored according to \code{ref}.
}
}
%\details{}
\value{
Returns a list with fields:
\item{x}{graph vertice coordinates.}
\item{p}{probability of each vertex belonging to each class, estimated from
MCMC samples in \code{cm}.}
}
\references{
For details on the specification of the degree-corrected stochastic
blockmodel, see "Bayesian Degree-Corrected Stochastic Blockmodels for
Community Detection".
}
%\author{}
%\note{}

\seealso{
\code{sbmlogit.map} and \code{sbmlogit.mcmc} to obtain \code{cm} lists.
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
