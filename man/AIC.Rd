\name{AIC}
\alias{AIC.randomLCA}
\title{AIC for randomLCA object}
\description{
Returns AIC for a randomLCA object.
}
\usage{
\method{AIC}{randomLCA}(object, \ldots, k = 2)
}
\arguments{
  \item{object}{randomLCA object}
  \item{\dots}{ additional argument; currently none is used. }
  \item{k}{penalty per parameter}
}
\value{
AIC
}
\examples{
\donttest{
pap.lca2 <- randomLCA(pap[, 1:10], freq = pap$Freq, nclass = 2)
print(AIC(pap.lca2))
}
}
\author{Ken Beath}
\keyword{methods}