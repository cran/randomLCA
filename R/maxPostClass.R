maxPostClass <- function(object)
  UseMethod("maxPostClass")

maxPostClass.randomLCA <- function(object){
  if (!inherits(object, "randomLCA"))
    stop("Use only with 'randomLCA' objects.\n")
  thepostprobs <- postClassProbs(object)
  maxprob <- apply(as.matrix(thepostprobs[,(dim(thepostprobs)[2]-object$nclass+1):dim(thepostprobs)[2]]),1,
                   function(x) max((x==max(x))*(1:object$nclas)))
  return(data.frame(thepostprobs[,1:(dim(thepostprobs)[2]-object$nclass)],maxProb=maxprob))
}
