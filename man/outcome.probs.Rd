\name{outcome.probs}
\alias{outcome.probs}
\alias{outcome.probs.randomLCA}
\title{Extract outcome probabilities for randomLCA object}
\description{
Extract outcome probabilities and confidence intervals for a randomLCA object.
}
\usage{
\method{outcome.probs}{randomLCA}(object, level = 0.95, boot=FALSE, type="perc",
	R=ifelse(type=="perc",999,200),\ldots)
}
\arguments{
  \item{object}{randomLCA object}
  \item{level}{confidence interval}
  \item{boot}{use parametric bootstrap to obtain confidence interval}
  \item{type}{type of bootstrap confidence intervals to use "perc", "normal" and "basic" are valid, see boot.ci for description. There is a risk of excessively large confidence limits when using other than "perc"}
  \item{R}{replications for parametric bootstrap}
  \item{\dots}{ additional argument; currently none is used. }
}
\details{
Confidence intervals are calculated based on asymptotic normality of the estimates transformed by either the inverse of the probit or logistic, or using parametric bootstrap. The asymptotic confidence intervals are currently only available for models without random effects.
}
\value{
}
\author{Ken Beath}
\examples{
# dentist data
data(dentistry)
# standard latent class with 2 classes
dentistry.lca2 <- randomLCA(dentistry[,1:5],freq=dentistry$freq,nclass=2)
print(outcome.probs(dentistry.lca2))
}
\keyword{methods}
