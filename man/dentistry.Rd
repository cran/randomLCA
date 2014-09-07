\name{dentistry}
\alias{dentistry}
\docType{data}
\title{Dental X-ray data}
\description{
  Six dentists evaluated dental x-rays for incipient caries in Handelman et al (1986), data consistss of 5 of the dentists analysed by Espeland and Handelman (1989) using a latent class model. Further analysis incorprating a random effects latent class model by Qu et al (1996), and by Albert and Dodd (2004)
}
\usage{data(dentistry)}
\format{
  A data frame with 32 observations on the following 6 variables.
  \describe{
    \item{\code{V1}}{Dentist 1}
    \item{\code{V2}}{Dentist 2}
    \item{\code{V3}}{Dentist 3}
    \item{\code{V4}}{Dentist 4}
    \item{\code{V5}}{Dentist 5}
    \item{\code{freq}}{Number of subjects}
  }
}
\source{
  Espeland and Handelman (1989)
  }
\references{
  Handelman, S.L., Leverett, D.H., Espeland, M.A. and Curzon, J.A. (1986) Clinical radiographic evaluation of sealed carious and sound tooth surfaces. \emph{Journal of the American Dental Association}, \bold{113}, 751--754.

  Espeland, M.A. and Handelman, S.L. (1989) Using latent class models to characterize and assess relative error in discrete distributions. \emph{Biometrics}, \bold{45}, 587--599.

  Qu, Y., Tan, M. and Kutner, M.H. (1996) Random effects models in latent
  class analysis for evaluating accuracy of diagnostic tests.
  \emph{Biometrics}, \bold{52}, 797--810.

  Albert P.S. and Dodd, L.E. (2004) A cautionary note on the robustness of Latent Class Models for estimating diagnostic error without a gold standard. \emph{Biometrics}, \bold{60}, 427--435.
}
\examples{
\dontrun{
data(dentistry)
# fit LCR model from Qu et al (1996)
dentistry.lca <- randomLCA(dentistry[,1:5],freq=dentistry$freq,level2size=5)
# start with constant laoding
dentistry.lcarandom <- randomLCA(dentistry[,1:5],freq=dentistry$freq,
	random=TRUE,probit=TRUE)
# allow loading to vary by dentist
dentistry.lcarandomunequal <- randomLCA(dentistry[,1:5],freq=dentistry$freq,
	random=TRUE,constload=FALSE,probit=TRUE)
}
}
\keyword{datasets}
