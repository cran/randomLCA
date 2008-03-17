\name{calc.cond.prob}
\alias{calc.cond.prob}
\title{Calculate Conditional Outcome Probabilities }
\description{
 Calculates the conditional outcome probabilities for random effects models or for standard latent class returns the outcome probabilities. For random effects, the outcome probabilities may be calculated for various levels of the random effect.
 }
\usage{
calc.cond.prob(object, conditionalp = 0.5)
}
\arguments{
  \item{object}{RandomLCA object}
  \item{conditionalp}{the p value for the random effect}
}
\value{
Returns a data frame containing class, block, outcome, outcomep (outcome probability) and perc (p value for random effect) if conditionalp is specified.
}
\author{Ken Beath \email{ken@kjbeath.com.au}
}
\keyword{methods}