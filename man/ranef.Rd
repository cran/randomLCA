\name{ranef}
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
A matrix with the first column containg the random effects and the second column the standard error of the random effects. 
}
\author{Ken Beath}
\examples{
}
\keyword{methods}