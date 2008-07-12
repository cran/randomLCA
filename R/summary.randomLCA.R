`summary.randomLCA` <-
function(object,...) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
	print(data.frame(Classes = object$nclass, AIC = object$aic, BIC = object$bic,
		log.Lik = c(object$log.Lik),row.names = " ") )
    cat("Class probabilities","\n")
	x <- object$classp
	names(x) <- paste("Class ",1:object$nclass)
	print(x,digits=4)
	if (object$random)  cat("Conditional outcome probabilities","\n")
    else cat("Outcome probabilities","\n")
	x <- as.data.frame(object$outcomep)
	row.names(x) <- paste("Class ",1:object$nclass)
	names(x) <- names(object$patterns)
	print(x,digits=4)
	if (object$random) {
		cat("Loadings","\n")
		x <- object$lambdacoef
		if (!object$byclass) x <- t(as.matrix(x))
		if (object$byclass) names1 <- paste("Class ",1:object$nclass)
		else names1 <- ''
		names2 <- 1:object$blocksize
		dimnames(x) <- list(names1,names2)
		print(x,digits=4)
		if (object$level2) {
			cat("Tau","\n")
			x <- object$taucoef
			if (object$byclass) {
				names(x) <- paste("Class ",1:object$nclass)
				print(x,digits=4)
			} else {
				x <- as.numeric(x)
				cat(x,"\n")
			}		
		}
	}
    out <- list()
	out$log.Lik <- object$log.Lik
	out$AIC <- object$aic
	out$BIC <- object$bic
	out$nclass <- object$nclass
	out$classp <- object$classp
	out$outcomep <- object$outcomep
	if (object$random) {
		out$lambdacoef <- object$lambdacoef
		if (object$level2) out$taucoef <- object$taucoef
	}
	class(out) <- "summ.randomLCA"
	invisible(out)
}

