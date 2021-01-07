\name{fitted}
\alias{fitted.randomLCA}
\title{fitted values}
\description{
Extract fitted values for randomLCA object.
}
\usage{
\method{fitted}{randomLCA}(object, \ldots)
}
\arguments{
  \item{object}{randomLCA object}
  \item{\dots}{ additional argument; currently none is used. }
}
\value{
A data frame. The first columns of the data frame correspond to the patterns, followed by the frequency of each pattern, and then the fitted number for each pattern.
}
\examples{
\donttest{
dentistry.lcarandom <- randomLCA(dentistry[, 1:5], freq = dentistry$freq, 
	random = TRUE, probit = TRUE, cores = 1)
print(fitted(dentistry.lcarandom))
}
}
\author{Ken Beath \email{ken.beath@mq.edu.au}
}
\keyword{methods}