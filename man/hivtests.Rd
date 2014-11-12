\name{hivtests}
\alias{hivtests}
\docType{data}
\title{HIV testing data}
\description{
  Serum samples are tested for HIV by 4 different biossays in Alvord et al (1988) and sensitivity and specificity determined using latent class analysis. Qu et al (1996) repeat the analysis using a model incorporating a random effect.
}
\usage{data(hivtests)}
\format{
  A data frame with 16 observations on the following 5 variables.
  \describe{
    \item{\code{V1}}{Test 1}
    \item{\code{V2}}{Test 2}
    \item{\code{V3}}{Test 3}
    \item{\code{V4}}{Test 4}
    \item{\code{freq}}{Number of subjects}
  }
}
\source{
  Qu, Tan and Kutner (1989)
  }
\references{
   Alvord, W.G., Drummond, J.E., Arthur, L.O., Goedert, J.J., Levine, P.H., Murphy, E.L., Weiss, S.H., and Blattner, W.A. (1988) A method for predicting individual HIV infection status in the absence of clinical information. \emph{AIDS Research and Human Retroviruses}, \bold{4}, 295--304.

  Qu, Y., Tan, M. and Kutner, M.H. (1996) Random effects models in latent
  class analysis for evaluating accuracy of diagnostic tests.
  \emph{Biometrics}, \bold{52}, 797--810.
}
\examples{
data(hivtests)
# fit standard latent class
hivtests.lca2 <- randomLCA(hivtests[,1:4],freq=hivtests$freq)
# with random effect model is not identifiable
}
\keyword{datasets}
