\name{calcMargProb}
\alias{calcMargProb}
\title{Calculates Marginal Outcome Probabilities}
\description{
Calculates the marginal outcome probabilities for a random effects latent class model, by integrating the outcome probability over the random effect. This is performed using Gauss-Hermite quadrature with the number of quadrature points specified for the model fitting.}
\usage{
calcMargProb(object)
}
\arguments{
  \item{object}{randomLCA object}
}
\value{
Returns a data frame containing class, block, outcome, outcomep (outcome probability).
}
\author{Ken Beath}
\keyword{methods}