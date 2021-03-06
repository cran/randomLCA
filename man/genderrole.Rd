\name{genderrole}
\alias{genderrole}
\docType{data}
\title{Gender Role Opinion Items}
\description{
 Opinions collected on gender roles in a study by Felling et al (1987). This was originally published in Heinen (1996) and subsequently in  Galindo Garre and Vermunt (2006).
}
\usage{genderrole}
\format{
  A data frame with 16 observations on the following 5 variables.
  \describe{
    \item{\code{Q1}}{Women's liberation sets women against men.}
    \item{\code{Q2}}{It's better for a wife not to have a job because that always poses problems in the household, especially if there are children.}
    \item{\code{Q3}}{The most natural situation occurs when the man is the breadwinner and the woman runs the household and takes care of the children.}
    \item{\code{Q4}}{It isn't really as important for a girl to get a good education as it is for a boy.}
    \item{\code{Q5}}{A woman is better suited to raise small children than a man.}
    \item{\code{Freq}}{Number of subjects}
  }
}
\source{
  Galindo Garre and Vermunt (2006)
  }
\references{
   Felling, A., Peters, J., and Schreuder, O. (1987) Religion in Dutch society 85: Documentation of a national survey on religious and secular attitudes in 1985. Amsterdam: Steinmetz Archive.
   
   Galindo Garre, F. and Vermunt, J.K. (2006) Avoiding boundary estimates in latent class analysis by Bayesian posterior mode estimation. \emph{Behaviormetrika}, \bold{33}, 43--59.
   
   Heinen, T. (1996) Latent Class and Discrete Latent Trait Models: Similarities and Differences. 
}
\examples{
\donttest{
# standard latent class
genderrole.lca1 <- randomLCA(genderrole[, 1:5], freq = genderrole$Freq, nclass = 1, cores = 1)
genderrole.lca2 <- randomLCA(genderrole[, 1:5], freq = genderrole$Freq, cores = 1)
genderrole.lca3 <- randomLCA(genderrole[, 1:5], freq = genderrole$Freq, nclass = 3, cores = 1)
# repeat with random effect with constant loading
# increase quadrature points and/or use higher penalty to obtain
# convergence
genderrole.lca1random <- randomLCA(genderrole[, 1:5], freq = genderrole$Freq, 
	nclass = 1, random = TRUE, cores = 1)
genderrole.lca2random <- randomLCA(genderrole[, 1:5], freq = genderrole$Freq, 
	random = TRUE, penalty = 0.1, quadpoints = 61, cores = 1)
genderrole.lca3random <- randomLCA(genderrole[, 1:5], freq = genderrole$Freq, 
	nclass = 3, random = TRUE, penalty = 0.1, quadpoints = 61, cores = 1)
# improved BIC for 1 class random
print(c(BIC(genderrole.lca1), BIC(genderrole.lca2), BIC(genderrole.lca3)))
print(c(BIC(genderrole.lca1random), BIC(genderrole.lca2random), 
	BIC(genderrole.lca3random)))
# can also repeat fits without constant loading to give mixture of IRT models
genderrole.lca1random2 <- randomLCA(genderrole[, 1:5], freq = genderrole$Freq, 
	nclass = 1, random = TRUE, constload = FALSE, cores = 1)
genderrole.lca2random2 <- randomLCA(genderrole[, 1:5], freq = genderrole$Freq, 
	random = TRUE, constload = FALSE, quadpoints = 61, penalty = 0.1, cores = 1)
genderrole.lca3random2 <- randomLCA(genderrole[, 1:5], freq = genderrole$Freq, 
	nclass = 3, random = TRUE, constload = FALSE, quadpoints = 61, penalty = 0.1, cores = 1)
# no improvement in fit
print(c(BIC(genderrole.lca1random2), BIC(genderrole.lca2random2), 
	BIC(genderrole.lca3random2)))
}
}
\keyword{datasets}
