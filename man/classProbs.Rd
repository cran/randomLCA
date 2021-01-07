\name{classProbs}
\alias{classProbs}
\title{
Determines class probabilities for fitted model
}
\description{
The class probabililities for the model are returned.
}
\usage{
classProbs(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{randomLCA object}
}
\details{
Simply extracts the corresponding variable from the randomLCA object.
}
\value{
A vector of class probabilities for each class.
}
\examples{
\donttest{
pap.lca2 <- randomLCA(pap[, 1:10], freq = pap$Freq, nclass = 2, cores = 1)
print(classProbs(pap.lca2))
}
}
\author{Ken Beath}
\keyword{methods}
