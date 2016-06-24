\name{randomLCA}
\alias{randomLCA}
\title{Fits a Latent Class Model including a Random Effect}
\description{
Fit latent class models, which may include a random effect.
}
\usage{
randomLCA(patterns,freq=NULL,nclass=2,calcSE=TRUE,notrials=20,
  random=FALSE,byclass=FALSE,quadpoints=21,constload=TRUE,blocksize=dim(patterns)[2],
  level2=FALSE,probit=FALSE,level2size=blocksize,
  qniterations=5,penalty=0.01,verbose=FALSE,
  seed = as.integer(runif(1, 0, .Machine$integer.max)))}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patterns}{Data frame or matrix of 0 and 1 defining the outcome patterns. May also include missing values, with randomLCA using maximum likelihood to fit the models using all available data.}
  \item{freq}{Frequency for each outcome pattern, if missing this is calculated from the patterns, and the patterns are summarised to remove duplicate values.}
  \item{nclass}{Number of classes to be fitted}
  \item{calcSE}{Calculate standard errors for parameters. Useful for bootstrapping.}
  \item{notrials}{For a standard latent class model, the number of random starting values used}
  \item{random}{Include random effect?}
  \item{byclass}{Vary random effect loading(s) by class?}
  \item{quadpoints}{Number of quadrature points for adaptive quadrature}
  \item{constload}{Outcome loadings are constant for random effects model?}
  \item{blocksize}{Where a random effects (single level) model is broken into blocks, that is the loadings are repeated, this defines the size of the blocks}
  \item{probit}{Probit model for random effect?}
  \item{level2}{Fit 2 level random effects model (further details to follow)?}
  \item{level2size}{Size of level 2 blocks if fitting 2 level models}
  \item{qniterations}{Number of Quasi-Newton iterations within each EM/adaptive cycle. Decrease if there is a failure to converge}
  \item{penalty}{penalty applied to likelihood for outcome probabilities. Shrinks outcome probabilities in slightly and can prevent extreme values. Setting penalty to 0 will produce an unpenalized fit.}
  \item{verbose}{Prints fit progress if true}
  \item{seed}{Initial random seed for generating starting values. This can be set to guarantee that the fit is the same each time, including the order of the classes.}
}
\details{
The structure of the patterns is assumed to be a number of blocks of different outcomes each of level2size, allowing outcomes to be repeated. Each outcome is assumed to have it's own loading.
An example is the width of the patterns is n and the level2size is n, resulting in n outcomes and therefore n loadings. Alternatively if the level2size is 1, then there are n repeats of the same outcome (but with different probabilities) with the same loading. In practice they may not be the same type of outcome, but usually will be.

The algorithm used is EM for the standard latent class and adaptive (in the sense of moving the location of the quadrature points) Gauss-Hermite quadrature for the random effects models. The number of quadrature points defaults to 21.
}
\value{
  randomLCA object
  This contains
\item{fit}{Fit object from optim}
\item{nclass}{Number of classes}
\item{classp}{Class probabilities}
\item{outcomep}{Outcome probability}
\item{lambdacoef}{Loadings}
\item{se}{Standard errors corresponding to results returned by optim}
\item{np}{Number of parameters}
\item{nobs}{Number of observations in total}
\item{logLik}{log likelihood for fitted model}
\item{penlogLik}{Penalised log likelihood for fitted model}
\item{observed}{Observed numbers corresponding to each pattern}
\item{fitted}{Fitted number corresponding to each pattern}
\item{deviance}{Deviance}
\item{classprob}{Posterior class probability for each pattern}
\item{bics}{BIC obtained for each trial when fitting initial latent class models}
\item{call}{call to randomLCA}
\item{random}{random parameter to randomLCA}
\item{constload}{constload parameter to randomLCA}
\item{level2}{level2 parameter to randomLCA}
\item{level2size}{level2size parameter to randomLCA}
\item{byclass}{byclass parameter to randomLCA}
\item{probit}{probit parameter to randomLCA}
\item{quadpoints}{quadpoints parameter to randomLCA}
\item{blocksize}{blocksize parameter to randomLCA}
\item{freq}{frequency of each pattern}
\item{qniterations}{qniterations parameter to randomLCA}
\item{penalty}{penalty parameter to randomLCA}
}
\note{In the returned object there a fields for patterns and frequencies. If frequencies are not supplied then the patterns and frequencies are constructed. If frequencies are supplied then zero rows are removed. When frequencies are supplied it is assumed that the data has been simplified. The returned fitted, posterior class probabilities etc, all correspond to the simplified patterns, not to the original data.

}
\author{Ken Beath}
\examples{
\dontrun{
# standard latent class with 2 classes
dentistry.lca2 <- randomLCA(dentistry[,1:5],freq=dentistry$freq,nclass=2)
# random effects model with constant random effect loading
dentistry.lca2random <- randomLCA(dentistry[,1:5],freq=dentistry$freq,
	nclass=2,random=TRUE,constload=TRUE,probit=TRUE)
# allow loading to vary by dentist
# this is the 2LCR model from Qu et al (1996)
dentistry.lca2random1 <- randomLCA(dentistry[,1:5],freq=dentistry$freq,
	nclass=2,random=TRUE,constload=FALSE,probit=TRUE)
}
}
\keyword{multivariate}
