`bestlca` <-
  function(patterns,freq,nclass,calcSE,notrials,probit,penalty,EMtol,verbose,cores) {
    
    
    thestatistic <- function(data, ...) {
      
      #browser()
      
      freq <- data[,1]
      patterns <- data[,2:dim(data)[2]]
      
      # noutliers <- max(1,round(dim(data)[1]*0.2))
      # outliers <- sample(c(rep(1,noutliers),rep(0,dim(data)[1]-noutliers)))
      
      #browser()
      thefit <- fitFixed(patterns,freq,nclass=nclass,initoutcomep=NULL,
               initclassp=NULL,calcSE=calcSE,justEM=FALSE,probit=probit,
               penalty=penalty,EMtol=EMtol,verbose=verbose, fullresults=FALSE)
      #browser()
    }
    
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    
    if (cores > 1) {
      if(.Platform$OS.type=="unix")  parallel <- "multicore"
      else parallel <- "snow" }
    else parallel <- "no"
        #browser()
    theboot <- boot(cbind(freq,patterns), thestatistic, R=notrials, sim = "parametric",
                    ran.gen = function(d, p) d,
                    parallel = parallel,
                    ncpus = cores,
                    calcSE=calcSE,justEM=TRUE,probit=probit,
                    penalty=penalty,EMtol=EMtol,verbose=verbose)
    #browser()
    res <- theboot$t
    if (verbose) for (i in 1:notrials) cat(c(res[i,2],res[[i]]$start.val),"\n")
    nfails <- sum(is.na(res[,2]))
    if (nfails > 0) warning(sprintf("Failed to obtain starting values for %i starting sets", nfails))
    
    bics <- -2*(res[,2])+log(res[,3])*res[,4]
    
    res <- res[order(res[,2], decreasing=TRUE),]
    res <- res[1,]
    
    
    #browser()

    classp <- res[5:(5+nclass-1)]
    outcomep <- matrix(res[(5+nclass):length(res)],nrow=nclass)
    
     maxlca <- fitFixed(patterns,freq,nclass=nclass,initoutcomep=outcomep,
                       initclassp=classp,calcSE=calcSE,justEM=FALSE,probit=probit,
                       penalty=penalty,EMtol=EMtol,verbose=verbose, fullresults=TRUE)
    if (verbose) {
      print("bic for class")
      print(bics)
    }
    return(c(maxlca,list(bics=bics)))
  }

