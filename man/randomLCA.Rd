\name{randomLCA}
\alias{randomLCA}
\title{Fits a Latent Class Model including a Random Effect}
\description{
Fit latent class models, which may include a random effect.
}
\usage{
randomLCA(patterns,freq,nclass=2,calcSE=FALSE,initmodel=NULL,blocksize=1,notrials=20,
	random=FALSE,byclass=FALSE,quadpoints=21,level2=FALSE,probit=FALSE,
	verbose=FALSE,seed = as.integer(runif(1, 0, .Machine$integer.max)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patterns}{Matrix of 0 and 1 definging the outcome patterns}
  \item{freq}{Frequency for each outcome pattern, if missing this is calculated from the patterns}
  \item{nclass}{Number of classes to be fitted}
  \item{calcSE}{Calculate standard errors for parameters}
  \item{initmodel}{Initial model of class randomLCA}
  \item{blocksize}{Size of blocks of outcomes}
  \item{notrials}{For a standard latent class model, the number of random starting values used}
  \item{random}{Random effect}
  \item{byclass}{Random effect by class}
  \item{quadpoints}{Number of quadrature points for adaptive quadrature}
  \item{level2}{Fit 2 level random effects model (further details to follow)}
  \item{probit}{Probit model for random effect}
  \item{seed}{Initial random seed for generating starting values}
  \item{verbose}{Prints verbose of fit progress}
}
\details{
The structure of the patterns is assumed to be a number of blocks of different outcomes each of blocksize, allowing outcomes to be repeated. Each outcome is assumed to have it's own loading.
An example is the width of the patterns is n and the blocksize is n, resulting in n outcomes and therefore n loadings. Alternatively if the blocksize is 1, then there are n repeats of the same outcome (but with different probabilities) with the same loading. In practice they may not be the same type of outcome, but usually will be.

The algorithm used is EM for the standard latent class and adaptive (in the sense of moving the location of teh quadrature points) Gauss-Hermite quadrature for the random effects models. The number of quadrature points defaults to 21.

There are some improvements that will be made to improve the speed. The first of these will be to allow Bayesian estimation or penalised likelihood using a Dirichlet prior for the outcome probabilities. This should allow a considerable improvement in convergence. The second will be to recode the major part of the likelihood calculations in C. Switching to a GEM algorithm to improve convergence may also be a possibility.
}
\value{
  randomLCA object
  This contains
\item{outcomep}{Outcome probability}
\item{classp}{Class probabilities}
\item{lambdacoef}{Loadings}
}
\author{Ken Beath}
\examples{
# dentist data
data(dentistry)
# standard latent class with 2 classes
dentistry.lca2 <- randomLCA(dentistry[,1:5],freq=dentistry$freq,nclass=2)
# random effects model with constant random effect loading
dentistry.lca2random <- randomLCA(dentistry[,1:5],freq=dentistry$freq,
	initmodel=dentistry.lca2,nclass=2,random=TRUE,probit=TRUE)
# allow loading to vary by dentist
# this is the 2LCR model from Qu et al (1996)
dentistry.lca2random1 <- randomLCA(dentistry[,1:5],freq=dentistry$freq,
	initmodel=dentistry.lca2random,nclass=2,random=TRUE,probit=TRUE,blocksize=5)
}
\keyword{multivariate}
