\name{plot}
\alias{plot.randomLCA}
\title{Plot a randomLCA object}
\description{
Plots the outcome probabilities for a randomLCA object, for random effects objects this can be either marginal or conditional or both. For a 2 level random effects model conditional2 will condition on the subject random effect and integrate over the period random effects. Note that plot is based on the xyplot function.
}
\usage{
\method{plot}{randomLCA}(x, \ldots ,graphtype = c("conditional", "marginal",
	"conditional2","both"), conditionalp = 0.5, classhorizontal = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{randomLCA object}
  \item{graphtype}{Type of graph}
  \item{conditionalp}{For a conditional graph the percentile corresponding to the random effect at which the outcome probability is to be calculated}
  \item{classhorizontal}{classes to be plotted across the page}
  \item{\dots}{additional parameters to xyplot}
}
\author{Ken Beath \email{ken@kjbeath.com.au}
}
\seealso{\code{\link{calc.cond.prob}},\code{\link{calc.cond.prob}}}
\examples{
data(uterinecarcinoma)
# standard latent class with 2 classes
uterinecarcinoma.lca2 <- randomLCA(uterinecarcinoma[,1:7],freq=uterinecarcinoma$freq)
plot(uterinecarcinoma.lca2)
}
\keyword{methods}