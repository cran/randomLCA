\name{calcCondProb}
\alias{calcCondProb}
\title{Calculate Conditional Outcome Probabilities }
\description{
 Calculates the conditional outcome probabilities for random effects models or for standard latent class returns the outcome probabilities. For random effects, the outcome probabilities may be calculated for various percentiles of the random effect.
 }
\usage{
calcCondProb(object, conditionalp = 0.5)
}
\arguments{
  \item{object}{RandomLCA object}
  \item{conditionalp}{the percentiles for the random effect}
}
\value{
Returns a data frame containing class, block, outcome, outcomep (outcome probability) and perc (percentiles of the random effect) if conditionalp is specified. For example a conditionalp of 0.5 is the 50th percentile or the median corresponding to a random effect of zero. 0.025 and 0.975 correspond to the 2.5th and 97.5th percential, so the region between them if 95\% of the variation in the data.
}
\author{Ken Beath \email{ken.beath@mq.edu.au}
}
\keyword{methods}