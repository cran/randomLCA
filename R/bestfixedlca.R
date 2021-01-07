`bestfixedlca` <-
  function(patterns,freq,nclass,calcSE,notrials,probit,penalty,EMtol,verbose,cores) {
    
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    
    bics <- rep(NA,notrials)
    maxll <- -Inf
    if (cores > 1) {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      res = foreach(i = 1:notrials, 
                    .options.RNG=seed[1]) %dorng% {
                      fitFixed(patterns,freq,nclass=nclass,calcSE=FALSE,justEM=TRUE,probit=probit,penalty=penalty,EMtol=EMtol,verbose=verbose)
                    }
      parallel::stopCluster(cl)
      for (i in 1:notrials) {
      if (verbose) cat(i," logLik = ",res[[i]]$logLik,"\n")
      bics[i] <- -2*(res[[i]]$logLik)+log(res[[i]]$nobs)*res[[i]]$np
      if (res[[i]]$logLik>maxll) {
        maxll <- res[[i]]$logLik
        maxlca <- res[[i]]
       }
    }
  } else {
    for (i in 1:notrials) {
      lca <- fitFixed(patterns,freq,nclass=nclass,calcSE=FALSE,justEM=TRUE,probit=probit,penalty=penalty,EMtol=EMtol,verbose=verbose)
       bics[i] <- -2*(lca$logLik)+log(lca$nobs)*lca$np
      #browser()
      if (lca$logLik > maxll) {
        maxll <- lca$logLik
        maxlca <- lca
      }
      if (verbose)
        cat("iteration ",i,"logLik = ",lca$logLik,"\n")
    }
  }
    if (verbose) print("refitting to obtain SE")
    #browser()
    maxlca <- fitFixed(patterns,freq,nclass=nclass,initoutcomep=maxlca$outcomep,
                       initclassp=maxlca$classp,calcSE=calcSE,justEM=FALSE,probit=probit,penalty=penalty,EMtol=EMtol,verbose=verbose)
    if (verbose) {
      print("bic for class")
      print(bics)
    }
    return(c(maxlca,list(bics=bics)))
  }

