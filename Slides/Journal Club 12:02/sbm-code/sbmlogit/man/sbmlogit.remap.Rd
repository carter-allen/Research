\name{sbmlogit.remap}
\alias{sbmlogit.remap}
\title{
  Remapping a configuration of labels to a canonical reference space.
}
\description{
  This routine remaps a configuration of labels \code{sigma} to a new vector
  \code{remap.sigma} such that if \code{ind(k)} is the index of the first
  appearance of label \code{k} in \code{remap.sigma} then
  \code{ind(i) < ind(j)} whenever \code{i < j}.
}
\usage{
sbmlogit.remap(sigma)
}
\arguments{
  \item{sigma}{
  a vector with numeric labels.
}
}
%\details{}

\value{
Returns a remapped \code{remap.sigma} for \code{sigma}.
}

\references{
For details on the specification of the degree-corrected stochastic
blockmodel, see "Bayesian Degree-Corrected Stochastic Blockmodels for
Community Detection".
}
%\author{}
%\note{}

%\seealso{}
\examples{
sigma <- c(2, 2, 3, 1, 3, 4, 2, 1)
remap.sigma <- sbmlogit.remap(sigma) # c(1, 1, 2, 3, 2, 4, 1, 3)
}

%\keyword{ ~kwd1 }
