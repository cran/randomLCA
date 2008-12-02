`bestlca` <-
function(patterns,freq,nclass,calcSE,notrials,probit,verbose,seed) {
	bics <- rep(NA,notrials)
    if(!exists(".Random.seed", envir = .GlobalEnv))
        runif(1)		     # initialize the RNG if necessary
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
	for (i in 1:notrials) {
		seed <- runif(1)*24674
		set.seed(seed)		
		lca <- fit.fixed.randomLCA(patterns,freq,nclass=nclass,calcSE=FALSE,probit=probit,verbose=verbose)
		currbic <- -2*(lca$logLik)+log(lca$nobs)*lca$np
		bics[i] <- currbic
		if (i==1) {
			maxbic <- currbic
			maxlca <- lca
		}
		if (currbic < maxbic) {
			maxbic <- currbic
			maxlca <- lca
		}
		if (verbose)
			cat("iteration ",i,"BIC ",BIC(lca),"\n")
	}
	if (calcSE) {
		if (verbose) print("refitting to obtain SE")
		maxlca <- fit.fixed.randomLCA(patterns,freq,nclass=nclass,initoutcomep=maxlca$outcomep,
			initclassp=maxlca$classp,calcSE=TRUE,probit=probit,verbose=verbose)
	}
	if (verbose) {
		print(maxlca)
		print("bic for class")
		print(bics)
	}
	assign(".Random.seed", RNGstate, envir = .GlobalEnv)
	return(c(maxlca,list(bics=bics)))
}

