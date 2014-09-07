\name{calc.cond2.prob}
\alias{calc.cond2.prob}
\title{Calculate Conditional Outcome Probabilities for 2 Level Models}
\description{
The conditional probabilities are obtained integrating over the period random effect.
 }
\usage{
calc.cond2.prob(object, conditionalp = 0.5)
}
\arguments{
  \item{object}{RandomLCA object}
  \item{conditionalp}{the p value for the subject random effect}
}
\value{
Returns a data frame containing perc (p value for random effect) , class, block, outcome and outcomep (outcome probability).
}
\author{Ken Beath \email{kenbeath@mq.edu.au}
}
\keyword{methods}