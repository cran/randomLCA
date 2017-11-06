\name{postClassProbs}
\alias{postClassProbs}
\title{
Determines posterior class probabilities for fitted model
}
\description{
The posterior class probabililities for each observed pattern and class is determined. These are returned as a data frame together with the patterns for each observation. If class=0 is requested then all classes are returned, otherwise only the selected class.
}
\usage{
postClassProbs(object,class=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{randomLCA object}
  \item{class}{class to be returned. Zero returns all classes.}
}
\details{
Extracts the corresponding data from the randomLCA object.
}
\value{
A data frame. The first columns of the data frame correspond to the patterns, followed by the frequency of each pattern, and then the posterior class probabilities for either the selected class or for all classes. The returned result is for the summarised data. If raw data is used, that is no frequencies, and it is required to calculated the posterior class probability for each observation then it is simply required to merge the class probabilities with the raw data, possibly removing any variable "freq" in the raw data. 
}
\examples{
pap.lca2 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=2)
print(postClassProbs(pap.lca2))
}
\author{Ken Beath}
\keyword{methods}
