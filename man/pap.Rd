\name{pap}
\alias{pap}
\docType{data}
\title{Positive Action program implementation}
\description{
The Positive Action program is a series of interventions designed to reduce negative behaviours in elementary-school students. In a study in Hawaii (Beets et al, 2006) information was recorded from students in the treatment group about whether the various parts of the program were implemented. While it is useful to describe the proportion of students experiencing implementation of each part of the program, latent class analysis will reveal if there are specific patterns to the implementation of the program (Alcock, 2008). 
}
\usage{pap}
\format{
  A data frame with 606 observations summarising the answers for 1566 students on the following 11 variables. For each variable the data has been dichotomized so that a 0 represents no implementation and a 1 represents some implementation. 
  \describe{
   \item{\code{Q1}}{you receive stickers from your teacher for doing positive actions?}
    \item{\code{Q2}}{you receive a word of the week card from your teacher?}
    \item{\code{Q3}}{you put notes in an icu box?}
    \item{\code{Q4}}{your teacher read notes about you from the icu box?}
    \item{\code{Q5}}{your teacher read your notes from the icu box?}
    \item{\code{Q6}}{your class receive a token for meeting your classroom goals?}
    \item{\code{Q7}}{you participate in a positive action assembly?}
    \item{\code{Q8}}{your class receive a balloon in an assembly for achieving their classroom goals?}
    \item{\code{Q9}}{your class participate in whole school positive action celebrations?}
    \item{\code{Q10}}{most weeks were you taught a positive action lesson?}
    \item{\code{Freq}}{Number of subjects}
  }
}
\source{
  Alcock (2008)
  }
\references{
Alcock, A. (2008). Latent Class Analysis (LCA) and Latent Profile Analysis (LPA). Retrieved from http://people.oregonstate.edu/~acock/growth/Mplus files/Latent Class Analysis Presentation.doc

Beets, M. W., Flay, B. R., Vuchinich, S., Snyder, F. J., Acock, A., Li, K. K.,  Durlak, J. (2009). Use of a social and character development program to prevent substance use, violent behaviors, and sexual activity among elementary-school students in Hawaii. American Journal of Public Health, 99(8), 1438-1445
}
\examples{
# standard latent class
pap.lca1 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=1)
pap.lca2 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=2)
\dontrun{
# standard latent class
pap.lca3 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=3)
pap.lca4 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=4)
pap.lca5 <- randomLCA(pap[,1:10],freq=pap$Freq,nclass=5)

# repeat with random effect with constant loading
# once BIC increases fitting further models is unnecessary
pap.lca1random <- randomLCA(pap[,1:10],
          freq=pap$Freq,nclass=1,random=TRUE)
pap.lca2random <- randomLCA(pap[,1:10],
          freq=pap$Freq,nclass=2,random=TRUE)
pap.lca3random <- randomLCA(pap[,1:10],
          freq=pap$Freq,nclass=3,random=TRUE)
# can also repeat fits without constant loading to give mixture of IRT models
pap.lca1random2 <- randomLCA(pap[,1:10],
          freq=pap$Freq,nclass=1,random=TRUE,constload=FALSE)

pap.lca2random2 <- randomLCA(pap[,1:10],
          freq=pap$Freq,nclass=2,random=TRUE,constload=FALSE)

pap.lca3random2 <- randomLCA(pap[,1:10],
          freq=pap$Freq,nclass=3,random=TRUE,constload=FALSE)

# produce table of BIC values
# shows 4 class best of standard latent class
# but 2 class latent class with constant loading has better BIC
pap.bic <- data.frame(bic=c(BIC(pap.lca1),BIC(pap.lca2),BIC(pap.lca3),BIC(pap.lca4),BIC(pap.lca5)),
                      bic2=c(BIC(pap.lca1random),BIC(pap.lca2random),BIC(pap.lca3random),NA,NA),
                      bic3=c(BIC(pap.lca1random2),BIC(pap.lca2random2),BIC(pap.lca3random2),NA,NA))
print(pap.bic)
# plot 4 class standard
plot(pap.lca4,type="b")
# plot 2 class standard
plot(pap.lca2random,type="b")
}
}
\keyword{datasets}
