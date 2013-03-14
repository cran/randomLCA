\name{uterinecarcinoma}
\alias{uterinecarcinoma}
\docType{data}
\title{Uterine Carcinoma Data}
\description{
 Classification of 118 histology samples by 118 pathologists. Original classification in Holmquist et al (1967) was to one of five categories, this has been reduced to two. Analysed by a number of authors, with a random effects in Qu et al (1996).
}
\usage{data(uterinecarcinoma)}
\format{
  A data frame with 20 observations on the following 8 variables.
  \describe{
    \item{\code{V1}}{Pathologist 1}
    \item{\code{V2}}{Pathologist 2}
    \item{\code{V3}}{Pathologist 3}
    \item{\code{V4}}{Pathologist 4}
    \item{\code{V5}}{Pathologist 5}
    \item{\code{V6}}{Pathologist 6}
    \item{\code{V7}}{Pathologist 7}
    \item{\code{freq}}{Number of observed pattern}
  }
}
\source{
  Qu et al (1996)
}
\references{

  Holmquist, N.D., McMahan, C.A., and Williams, O.D. (1967) Variability in classificationof carcinoma in situ of the uterine cervix.
  \emph{Archives of Pathology}, \bold{84}, 344--345.
  Qu, Y., Tan, M. and Kutner, M.H. (1996) Random effects models in latent
  class analysis for evaluating accuracy of diagnostic tests.
  \emph{Biometrics}, \bold{52}, 797--810.
}
\examples{
\dontrun{
data(uterinecarcinoma)
uterinecarcinoma.lca2 <- randomLCA(uterinecarcinoma[,1:7],
	freq=uterinecarcinoma$freq)

uterinecarcinoma.lcarandom2 <- randomLCA(uterinecarcinoma[,1:7],
	freq=uterinecarcinoma$freq,initmodel=uterinecarcinoma.lca2,
	random=TRUE,probit=TRUE,quadpoints=41)
# LCR1 model of Que et al. This is fairly unstable and requires
#	starting values from the simpler model without loadings by class
# it is also slow and doesn't improve the model fit
uterinecarcinoma.lcarandom2by <- randomLCA(uterinecarcinoma[,1:7],freq=uterinecarcinoma$freq,
	initmodel=uterinecarcinoma.lcarandom2,byclass=TRUE,
	random=TRUE,probit=TRUE,quadpoints=61)}
}
\keyword{datasets}
