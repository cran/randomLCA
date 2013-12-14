\name{genderrole}
\alias{genderrole}
\docType{data}
\title{Gender Role Opinion Items}
\description{
 Opinions collected on gender roles in a study by Felling et al (1987). This was originally published in Heinen (1996) and subsequently in  Galindo Garre and Vermunt (2006).
}
\usage{data(genderrole)}
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
\dontrun{
data(genderrole)
# standard latent class
genderrole.lca1 <- randomLCA(genderrole[,1:5],freq=genderrole$Freq,nclass=1)
genderrole.lca2 <- randomLCA(genderrole[,1:5],freq=genderrole$Freq)
genderrole.lca3 <- randomLCA(genderrole[,1:5],freq=genderrole$Freq,nclass=3)
# repeat with random effect
genderrole.lca1random <- randomLCA(genderrole[,1:5],freq=genderrole$Freq,
	nclass=1,initmodel=genderrole.lca1,random=TRUE)
genderrole.lca2random <- randomLCA(genderrole[,1:5],freq=genderrole$Freq,
	random=TRUE,initmodel=genderrole.lca2,quadpoints=41)
genderrole.lca3random <- randomLCA(genderrole[,1:5],freq=genderrole$Freq,
	nclass=3,random=TRUE,initmodel=genderrole.lca3,quadpoints=51)
# improved BIC for 1 class random
print(c(BIC(genderrole.lca1),BIC(genderrole.lca2),BIC(genderrole.lca3)))
print(c(BIC(genderrole.lca1random),BIC(genderrole.lca2random),
	BIC(genderrole.lca3random)))
# can also repeat fits with blocksize=5 to give mixture of IRT models
genderrole.lca1random2 <- randomLCA(genderrole[,1:5],freq=genderrole$Freq,
	nclass=1,initmodel=genderrole.lca1,random=TRUE,blocksize=5)
genderrole.lca2random2 <- randomLCA(genderrole[,1:5],freq=genderrole$Freq,
	random=TRUE,initmodel=genderrole.lca2,blocksize=5,quadpoints=41)
genderrole.lca3random2 <- randomLCA(genderrole[,1:5],freq=genderrole$Freq,
	nclass=3,random=TRUE,initmodel=genderrole.lca3,blocksize=5,quadpoints=51)
	# no improvement in fit
print(c(BIC(genderrole.lca1random2),BIC(genderrole.lca2random2),
	BIC(genderrole.lca3random2)))
}
}
\keyword{datasets}
