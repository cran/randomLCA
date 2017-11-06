\name{logLik}
\alias{logLik.randomLCA}
\title{log Likelikelihood for randomLCA object}
\description{
Returns log Likelihood for a randomLCA object.
}
\usage{
\method{logLik}{randomLCA}(object, \ldots)
}
\arguments{
  \item{object}{randomLCA object}
  \item{\dots}{ additional argument; currently none is used. }
}
\value{
The loglikelihood.
}
\examples{
pap.lca2 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=2)
print(logLik(pap.lca2))
}
\author{Ken Beath}
\keyword{methods}