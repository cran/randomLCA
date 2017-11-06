\name{simulate}
\alias{simulate.randomLCA}
\title{Simulate}
\description{
Simulate data from a fitted randomLCA model}
\usage{
\method{simulate}{randomLCA}(object, nsim, seed, \ldots)
}
\arguments{
  \item{object}{randomLCA object}
  \item{nsim}{number of data sets to be simulated}
  \item{seed}{random seed}
  \item{\ldots}{additional optional arguments.}
}
\details{
Generates random data from the supplied object.
}
\value{
A simulated data frame or a list of simulated data frames.
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
print((sum(simlrt>=obslrt)+1)/(nsims+1))}
}
\author{Ken Beath}
\keyword{methods}
