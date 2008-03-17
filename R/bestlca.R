`bestlca` <-
function(patterns,freq,nclass,calcSE,notrials,verbose) {
	bics <- rep(NA,notrials)
	set.seed(24567)
	for (i in 1:notrials) {
		seed <- runif(1)*24674
		set.seed(seed)		
		lca <- fit.fixed.randomLCA(patterns,freq,nclass=nclass,calcSE=FALSE,verbose=verbose)
		bics[i] <- lca$bic
		if (i==1) {
			maxbic <- lca$bic
			maxlca <- lca
		}
		if (lca$bic < maxbic) {
			maxbic <- lca$bic
			maxlca <- lca
		}
		if (verbose)
			cat("iteration ",i,"BIC ",lca$bic,"\n")
	}
	if (calcSE) {
		if (verbose) print("refitting to obtain SE")
		maxlca <- fit.fixed.randomLCA(patterns,freq,nclass=nclass,initoutcomep=maxlca$outcomep,
			initclassp=maxlca$classp,calcSE=TRUE,verbose=verbose)
	}
	if (verbose) {
		print(maxlca)
		print("bic for class")
		print(bics)
	}
	return(c(maxlca,list(bics=bics)))
}

