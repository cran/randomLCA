\name{ranef}
\alias{ranef}
\alias{ranef.randomLCA}
\title{Extract random effects from a randomLCA object}
\description{
Extracts the Empirical Bayes estimates of the random effects.
}
\usage{
\method{ranef}{randomLCA}(object, \ldots)
}
\arguments{
  \item{object}{randomLCA object with a random effect}
  \item{\dots}{ additional argument; currently none is used. }
}
\value{
A matrix with the first column containing the random effects and the second column the standard error of the random effects. 
}
\examples{
\donttest{
pap.lca2random <- randomLCA(pap[, 1:10], freq = pap$Freq, random = TRUE, nclass = 2, cores = 1)
print(ranef(pap.lca2random))
}
}
\author{Ken Beath}
\keyword{methods}