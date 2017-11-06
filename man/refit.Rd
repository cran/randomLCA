\name{refit}
\alias{refit}
\alias{refit.randomLCA}
\title{Refit an randomLCA object}
\description{
Refits an randomLCA object using new data. This is useful when fitting simulated data, for example using a bootstrap. 
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
The useunit parameter determines whether the parameter estimates from the supplied model are used as initial values or whether the complete model fitting process is repeated. If the initial values are used then fitting will be faster, and the fitted classes will be similar to those in the original model. If the data was not generated from the original model there is an increased risk that the fit will not find the global maxima. For this reason when performing a bootstrap Likelihood ratio test it is better to use useinit=FALSE. However when using useinit=FALSE there may be label switching, where the estimated classes are similar, but in a different order. Unless the estimated parameters are assigned to the correct classes this will invalidate the results of a parametric bootstrap for parameter confidence intervals.
}
\examples{
\dontrun{
myocardial.lca1 <- randomLCA(myocardial[,1:4],freq=myocardial$freq,nclass=1)
myocardial.lca2 <- randomLCA(myocardial[,1:4],freq=myocardial$freq)
# calculate observed lrt
obslrt <- 2*(logLik(myocardial.lca2)-logLik(myocardial.lca1))

print(obslrt)

nsims <- 999
# generate the simulations
thesims <- simulate(myocardial.lca1, nsims)
# for each simulation determin lrt
simlrt <- rep(NA,nsims)
for (isim in 1:nsims) {
  submodel <- refit(myocardial.lca1,newpatterns=thesims[[isim]])
  fullmodel <- refit(myocardial.lca2,newpatterns=thesims[[isim]])
  simlrt[isim] <- 2*(logLik(fullmodel)-logLik(submodel))
  print(c(isim,simlrt[isim]))
}
# calculate p value as proportion of simulated lrt greater than observed,
# corrected by adding one to numerator and denominator
print((sum(simlrt>=obslrt)+1)/(nsims+1))
}
}
\author{Ken Beath}
\keyword{methods}