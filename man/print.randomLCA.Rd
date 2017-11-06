\name{print.randomLCA}
\alias{print.randomLCA}
\title{print for randomLCA object}
\description{
Prints a randomLCA object. Prints summary.
}
\usage{
\method{print}{randomLCA}(x, \ldots)
}
\arguments{
  \item{x}{randomLCA object}
  \item{\dots}{ additional argument; currently none is used. }
}
\examples{
pap.lca1 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=1)
pap.lca1
# or
print(pap.lca1)
}
\author{Ken Beath}
\keyword{methods}