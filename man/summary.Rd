\name{summary}
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
\details{
Not complete.
}
\value{
  \item{logLik}{Log Likelihood}
  \item{AIC}{AIC}
  \item{BIC}{BIC}
  \item{nclass}{no of classes}
  \item{probit}{link is probit}
  \item{classp}{class probabilities}
  \item{outcomep}{outcome probabilities (conditional)}
  \item{margoutcomep}{outcome probabilities (marginal), if model contains random effects}
  \item{random}{model includes random effects}
  \item{level2}{model has 2 level hierarchy}
  \item{byclass}{lambda and tau vary by class}
  \item{lambdacoef}{lambda coefficients}
  \item{taucoef}{tau coefficients}
}
\author{Ken Beath}
\examples{
}
\keyword{methods}