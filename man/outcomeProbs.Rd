\name{outcomeProbs}
\alias{outcomeProbs}
\alias{outcomeProbs.randomLCA}
\title{Extract outcome probabilities for randomLCA object}
\description{
Extract outcome probabilities and confidence intervals for a randomLCA object.
}
\usage{
\method{outcomeProbs}{randomLCA}(object,level = 0.95, boot=FALSE, type="norm",
    R=ifelse(type=="norm",199,999),scale=c("prob","raw"),\ldots)
}
\arguments{
  \item{object}{randomLCA object}
  \item{level}{confidence interval}
  \item{boot}{use parametric bootstrap to obtain confidence interval}
  \item{type}{type of bootstrap confidence intervals to use, with "perc" or "norm" valid, see boot.ci for description. It seems reasonable to use the normal approximation.}
  \item{R}{replications for parametric bootstrap}
  \item{scale}{either "prob" where probabilities are returned, the default, or "raw" where the probabilities are returned on the logit or probit scale, depending on which scale was selected in the texttt{randomLCA} function}
  \item{\dots}{ additional argument; currently none is used. }
}
\details{
Confidence intervals are calculated based on asymptotic normality of the estimates transformed by either the inverse of the probit or logistic, or using parametric bootstrap. The asymptotic confidence intervals are currently only available for models without random effects. For the confidence intervals obtained from the parametric bootstrap, the bootstrap is performed on the data that has been transformed to the logit or probit scale, as appropriate. The samples are close to normal allowing for the use of confidence intervals based on the normal approximation. About 199 replications gives similar accuracy to percentile with 999.
}
\value{
Data frame consisting of outcome probabilities and confidence intervals. One for each class.
}
\author{Ken Beath}
\examples{
# standard latent class with 2 classes
dentistry.lca2 <- randomLCA(dentistry[,1:5], freq=dentistry$freq, nclass=2)
print(outcomeProbs(dentistry.lca2))
# print on the default logit scale
print(outcomeProbs(dentistry.lca2, scale="raw"))
# convert back to probabilities
print(1.0/(1.0+exp(-outcomeProbs(dentistry.lca2, scale="raw")[[1]])))
print(1.0/(1.0+exp(-outcomeProbs(dentistry.lca2, scale="raw")[[2]])))
}
\keyword{methods}
