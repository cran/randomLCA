`bestlca` <-
function(patterns,freq,nclass,calcSE,notrials,probit,verbose) {
	bics <- rep(NA,notrials)
	for (i in 1:notrials) {
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
	return(c(maxlca,list(bics=bics)))
}

