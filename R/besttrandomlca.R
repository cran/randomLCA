bestrandomlca <-
  function(patterns,freq,nclass,calcSE,notrials,probit,constload,byclass,blocksize,quadpoints,penalty,qniterations,EMtol,verbose,cores) {
    
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    
    bics <- rep(NA,notrials)
    maxll <- -Inf
    nfail <- 0
    if (cores > 1) {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      res = foreach(i = 1:notrials, 
                    .packages = c("fastGHQuad","Rcpp"),
                     .options.RNG=seed[1]) %dorng% {
                      fitAdaptRandom(patterns,freq=freq,nclass=nclass,calcSE=FALSE,justEM=TRUE,
                                     initoutcomep=NULL,initclassp=NULL,initlambdacoef=NULL,
                                     initmomentdata=NULL,
                                     gh=norm.gauss.hermite(quadpoints),
                                     constload=constload,probit=probit,
                                     byclass=byclass,blocksize=blocksize,
                                     qniterations=qniterations,
                                     penalty=penalty,EMtol=EMtol,verbose=FALSE)
                    }
      parallel::stopCluster(cl)
      for (i in 1:notrials) {
        if (is.null(res[[i]])) nfail <- nfail+1
        else {
          if (verbose) cat(i," logLik = ",res[[i]]$logLik,"\n")
          bics[i] <- -2*(res[[i]]$logLik)+log(res[[i]]$nobs)*res[[i]]$np
          if (res[[i]]$logLik>maxll) {
            maxll <- res[[i]]$logLik
            maxlca <- res[[i]]
          }
        }
      }
    } else {
      for (i in 1:notrials) {
        lca <-   fitAdaptRandom(patterns,freq=freq,nclass=nclass,calcSE=FALSE,justEM=TRUE,
                                initoutcomep=NULL,initclassp=NULL,initlambdacoef=NULL,
                                initmomentdata=NULL,
                                gh=norm.gauss.hermite(quadpoints),
                                constload=constload,probit=probit,
                                byclass=byclass,blocksize=blocksize,
                                qniterations=qniterations,
                                penalty=penalty,EMtol=EMtol,verbose=verbose)
        if (is.null(lca)) nfail <- nfail+1
        else {
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
    }
     if (nfail > 0) warning(sprintf("Covengence failed for %i starting values. Try increasing quadrature points.",nfail))
    if (maxll==-Inf) stop("No starting values found. Try increasing quadrature points.")
      if (verbose) print("refitting")
      maxlca <- fitAdaptRandom(patterns,freq=freq,nclass=nclass,calcSE=calcSE,justEM=FALSE,
                               initoutcomep=maxlca$outcomep,initclassp=maxlca$classp,
                               initlambdacoef=maxlca$lambdacoef,
                               initmomentdata=maxlca$momentdata,
                               gh=norm.gauss.hermite(quadpoints),
                               constload=constload,probit=probit,
                               byclass=byclass,blocksize=blocksize,
                               qniterations=qniterations,
                               penalty=penalty,EMtol=EMtol,verbose=verbose)
     # browser()
    return(c(maxlca,list(bics=bics)))
  }

