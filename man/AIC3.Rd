\name{AIC3}
\alias{AIC3}
\title{AIC with 3 penalty for randomLCA object}
\description{
Returns AIC with penalty 3 for a randomLCA object.
}
\usage{
AIC3(object)
}
\arguments{
  \item{object}{randomLCA object}
}
\value{
AIC3.
}
\examples{
\donttest{
pap.lca2 <- randomLCA(pap[, 1:10], freq = pap$Freq, nclass = 2, cores = 1)
print(AIC3(pap.lca2))
}
}
\author{Ken Beath}
\keyword{methods}