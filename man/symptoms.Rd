\name{symptoms}
\alias{symptoms}
\docType{data}
\title{Symptoms data}
\description{
  Allergy and respiratory symptoms for infants 0 to 2 years in six month periods. Outcome is presence or absence of symptom in the six months. Original data was collected at Visits 1-7 over the 2 year period which were summarised to six month periods.
  
  This is data relating to a paper which will hopefully be published shortly. I will be updating this entry then.
  
  Note that these models can be slow to fit, with the "symptoms.lca2random2" model taking about 3 hours. I have a couple more possibilities to speed things up. These are to code the algorithms completely in C (they are currently a mix of C and R), use a Bayesian methodology to prevent the estimates moving towards the boundaries (this can also be described as penalized likelihood) and to look at other optimisation functions.
  
  Thanks to the investigators of the CAPS study for making the data available.
}
\usage{data(symptoms)}
\format{
  A data frame with 444 observations on the following 17 variables.
  \describe{
	\item{\code{Nightcough.13}}{Night cough in visits 1-3}
	\item{\code{Wheeze.13}}{Wheeze in visits 1-3}
	\item{\code{Itchyrash.13}}{Itchy rash in visits 1-3}
	\item{\code{FlexDerma.13}}{Flexural Dermatitis in visits 1-3}
	\item{\code{Nightcough.45}}{Night cough in visits 1-3}
	\item{\code{Wheeze.45}}{Wheeze in visits 4-5}
	\item{\code{Itchyrash.45}}{Itchy rash in visits 4-5}
	\item{\code{FlexDerma.45}}{Flexural Dermatitis in visits 4-5}
	\item{\code{Nightcough.6}}{Night cough in visit 6}
	\item{\code{Wheeze.6}}{Wheeze in visit 6}
	\item{\code{Itchyrash.6}}{Itchy rash in visit 6}
	\item{\code{FlexDerma.6}}{Flexural Dermatitis in visits 1-3}
	\item{\code{Nightcough.7}}{Night cough in visit 7}
	\item{\code{Wheeze.7}}{Wheeze in visit 7}
	\item{\code{Itchyrash.7}}{Itchy rash in visit 7}
	\item{\code{FlexDerma.7}}{Flexural Dermatitis in visit 7}
	\item{\code{Freq}}{Number of subjects}
   }
}
\source{
  Mihrshai et al (2001)
  }
\references{
   Mihrshahi, S.,  Peat, J.K., Webb, K., Tovey, R.E., Marks, G.B., Mellis, C.M. and Leeder S.R. (2001) The Childhood Asthma Prevention Study (CAPS): Design and research protocol of a randomized trial for the primary prevention of asthma. Control led Clinical Trials, \bold{22}:333--354.

}
\examples{
data(symptoms)
\dontrun{
symptoms.lca2 <- randomLCA(symptoms[,1:16],freq=symptoms$Freq,nclass=2)
symptoms.lca2random <- randomLCA(symptoms[,1:16],freq=symptoms$Freq,initmodel=symptoms.lca2,random=TRUE,nclass=2,blocksize=4)
symptoms.lca2random2 <- randomLCA(symptoms[,1:16],freq=symptoms$Freq,initmodel=symptoms.lca2random,random=TRUE,level2=TRUE,nclass=2,blocksize=4)
}
}
\keyword{datasets}
