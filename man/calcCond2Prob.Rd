\name{calcCond2Prob}
\alias{calcCond2Prob}
\title{Calculate Conditional Outcome Probabilities for 2 Level Models}
\description{
The conditional probabilities are obtained integrating over the period random effect.
 }
\usage{
calcCond2Prob(object, conditionalp = 0.5)
}
\arguments{
  \item{object}{RandomLCA object}
  \item{conditionalp}{the percentiles for the random effect}
}
\value{
Returns a data frame containing class, block, outcome, outcomep (outcome probability) and perc (percentiles of the random effect) if conditionalp is specified.  For example a conditionalp of 0.5 is the 50th percentile or the median corresponding to a random effect of zero. 0.025 and 0.975 correspond to the 2.5th and 97.5th percential, so the region between them if 95\% of the variation in the data.
}
\author{Ken Beath \email{kenbeath@mq.edu.au}
}
\examples{
\dontrun{
symptoms.lca2random2 <- randomLCA(symptoms[,1:16],freq=symptoms$Freq,
	random=TRUE,level2=TRUE,nclass=2,level2size=4,constload=FALSE)
print(calcCond2Prob(symptoms.lca2random2))
}
}
\keyword{methods}