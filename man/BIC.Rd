\name{BIC}
\alias{BIC.randomLCA}
\title{BIC for randomLCA object}
\description{
Returns BIC for a randomLCA object.
}
\usage{
\method{BIC}{randomLCA}(object, \ldots)
}
\arguments{
  \item{object}{randomLCA object}
  \item{\dots}{ additional argument; currently none is used. }
}
\value{
BIC
}
\examples{
pap.lca2 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=2)
print(BIC(pap.lca2))
}
\author{Ken Beath}
\keyword{methods}