\name{outcome.probs}
\alias{outcome.probs}
\alias{outcome.probs.randomLCA}
\title{Extract outcome probabilities for randomLCA object}
\description{
Extract outcome probabilities and confidence intervals for a randomLCA object.
}
\usage{
\method{outcome.probs}{randomLCA}(object, level = 0.95, \ldots)
}
\arguments{
  \item{object}{randomLCA object}
  \item{level}{confidence interval}
  \item{\dots}{ additional argument; currently none is used. }
}
\details{
Confidence intervals are calculated based on asymptotic normality of the estimates transformed by either the inverse of the probit or logistic.
}
\value{
}
\author{Ken Beath}
\examples{
}
\keyword{methods}
