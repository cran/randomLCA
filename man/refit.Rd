\name{refit}
\alias{refit}
\alias{refit.randomLCA}
\title{Refit an randomLCA object}
\description{
Refits an randomLCA object using new data. For an example, see the simulate method.
}
\usage{
\method{refit}{randomLCA}(object, newpatterns, newfreq, useinit=FALSE, \ldots)
}
\arguments{
  \item{object}{randomLCA object}
  \item{newpatterns}{the new patterns that are to be fitted using the existing model}
  \item{newfreq}{the frequencies corresponding to the patterns if required}
  \item{useinit}{use initial values from randomLCA object}
  \item{\dots}{ additional argument; currently none is used. }
}
\value{
The fitted model to the new data.
}
\details{
The useunit parameter determines whether the parameter estimates from the supplied model are used as initial values or whether the complete model fitting process is repeated. If the initial values are used then fitting will be faster, and the fitted classes will be similar to those in the original model. If the data was not generated from the original model there is an increased risk that the fit will not find the global maxima. For this reason when performing a bootstrap Likelihood ratio test it is better to use useinit=FALSE. However when using useinit=FALSE there may be label switching, where the estimated classes are similar, but in a different order. Unless the estimated parameters are assigned to teh correct classes this will invalidate the results of a parametric bootstrap for parameter confidence intervals.
}
\author{Ken Beath}
\keyword{methods}