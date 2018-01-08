\name{maxPostClass}
\alias{maxPostClass}
\alias{maxPostClass.randomLCA}
\title{
Determines class with maximum posterior class probability for each observation
}
\description{
For each observation the posterior class probability is determined for each class, and then the class with the maximum posterior class probability is returned.
}
\usage{
maxPostClass(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{randomLCA object}
}
\details{
Returns the class with the maximum posterior class probability for each observation.
}
\value{
A data frame. The first columns of the data frame correspond to the patterns, followed by the frequency of each pattern, and then the class with the maximum posterior class probability. The returned result is for the summarised data. If raw data is used, that is no frequencies, and it is required to calculated the posterior class probability for each observation then it is simply required to merge the maximum class with the raw data, possibly removing any variable "Freq" in the raw data. 
}
\examples{
pap.lca2 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=2)
maxClass.lca2 <- maxPostClass(pap.lca2)
names(maxClass.lca2)[length(names(maxClass.lca2))] <- "maxProb.lca2"

pap.lca3 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=3)
maxClass.lca3 <- maxPostClass(pap.lca3)
names(maxClass.lca3)[length(names(maxClass.lca3))] <- "maxProb.lca3"

maxClass <- merge(maxClass.lca2,maxClass.lca3)
# aggregate because there is more than one record for each combination
# shows the relationship between the 2 and 3 class models
print(aggregate(maxClass$Freq,list(maxProb2=maxClass$maxProb.lca2,
    maxProb3=maxClass$maxProb.lca3),sum))
}
\author{Ken Beath}
\keyword{methods}
