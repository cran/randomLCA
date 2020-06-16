\name{summary.randomLCA}
\alias{summary.randomLCA}
\title{Summary for randomLCA object}
\description{
Summarises the fit of a randomLCA object.
}
\usage{
\method{summary}{randomLCA}(object, \ldots)
}
\arguments{
  \item{object}{randomLCA object}
  \item{\dots}{ additional argument; currently none is used. }
}
\value{
  \item{logLik}{Log Likelihood}
  \item{AIC}{AIC}
  \item{BIC}{BIC}
  \item{AIC3}{AIC with penalty of 3}
  \item{nclass}{no of classes}
  \item{probit}{link is probit}
  \item{classp}{class probabilities}
  \item{outcomep}{outcome probabilities (conditional)}
  \item{margoutcomep}{outcome probabilities (marginal), if model contains random effects}
  \item{random}{model includes random effects}
  \item{level2}{model has 2 level hierarchy}
  \item{constload}{loadings are constant by outcome}
  \item{byclass}{lambda and tau vary by class}
  \item{lambdacoef}{lambda coefficients}
  \item{taucoef}{tau coefficients}
}
\examples{
\donttest{
pap.lca1 <- randomLCA(pap[, 1:10], freq = pap$Freq, nclass = 1)
summary(pap.lca1)
}
}
\author{Ken Beath}
\keyword{methods}