\name{plot}
\alias{plot.randomLCA}
\title{Plot a randomLCA object}
\description{
Plots the outcome probabilities for a randomLCA object, for random effects objects this can be either marginal or conditional or both. For a 2 level random effects model conditional2 will condition on the subject random effect and integrate over the period random effects. Note that plot is based on the xyplot function.
}
\usage{
\method{plot}{randomLCA}(x, \ldots , graphtype=ifelse(x$random, "marginal","conditional"),
     conditionalp=0.5,classhorizontal=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{randomLCA object}
  \item{graphtype}{Type of graph}
  \item{conditionalp}{For a conditional graph the percentile corresponding to the random effect at which the outcome probability is to be calculated}
  \item{classhorizontal}{classes to be plotted across the page}
  \item{\dots}{additional parameters to xyplot}
}
\author{Ken Beath \email{ken.beath@mq.edu.au}
}
\seealso{\code{\link{calcCondProb}},\code{\link{calcMargProb}}}
\examples{
\dontrun{
# standard latent class with 2 classes
uterinecarcinoma.lca2 <- randomLCA(uterinecarcinoma[,1:7],freq=uterinecarcinoma$freq)
plot(uterinecarcinoma.lca2)
uterinecarcinoma.lcarandom2 <- randomLCA(uterinecarcinoma[,1:7],
  freq=uterinecarcinoma$freq,random=TRUE,probit=TRUE,quadpoints=41)
# default for random effects models is marginal
plot(uterinecarcinoma.lcarandom2)
# default for random effects models conditional is p=0.5 i.e. median
plot(uterinecarcinoma.lcarandom2,graphtype="conditional")
# look at variability by plotting conditional probabilities at 0.05,0.5 and 0.95
plot(uterinecarcinoma.lcarandom2,graphtype="conditional",conditionalp=c(0.05,0.5,0.95))
}
}
\keyword{methods}