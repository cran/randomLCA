\name{myocardial}
\alias{myocardial}
\docType{data}
\title{Myocardial Infarction (MI)}
\description{
	Four tests were performed on hospital patients to determine if a myocardial infarction had occurred. 
  }
\usage{myocardial}
\format{
  A data frame with 32 observations on the following 6 variables.
  \describe{
    \item{\code{Q.wave}}{result from ECG test}
    \item{\code{History}}{clinical history}
    \item{\code{LDH}}{flipped, enzyme related to tissue breakdown}
    \item{\code{CPK}}{high, creatine kinase or creatine phosphokinase, related to muscle damage}
    \item{\code{freq}}{Number of subjects}
  }
}
\source{
  Rindskopf and Rindskopf (1986)
  }
\references{
  Galen, R.S. and Gambino, S.R. (1975) Beyond Normality: The Predictive Value and Efficiency of Medical Diagnosis. Wiley:New York.

  Rindskopf, D. and Rindskopf, W. (1986) The Value of Latent Class Analysis in Medical Diagnosis. \emph{Statistics in Medicine}, \bold{5}, 21--27.
}
\examples{
\donttest{
# fit 2 class model from Rindskopf and Rindskopf (1986)
myocardial.lca2 <- randomLCA(myocardial[, 1:4], freq = myocardial$freq, cores = 1)
}
}
\keyword{datasets}
