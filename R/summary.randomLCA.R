summary.randomLCA <-
function(object,...) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
    out <- list()
	out$log.Lik <- object$log.Lik
	out$AIC <- object$aic
	out$BIC <- object$bic
	out$nclass <- object$nclass
	out$classp <- object$classp
	names(out$classp) <- paste("Class ",1:object$nclass)
	out$outcomep <- as.data.frame(object$outcomep)
	row.names(out$outcomep) <- paste("Class ",1:object$nclass)
	names(out$outcomep) <- names(object$patterns)
	out$random <- object$random
	out$level2 <- object$level2
	if (object$random) {
		out$lambdacoef <- object$lambdacoef
		if (!object$byclass) out$lambdacoef <- t(out$lambdacoef)
		if (object$byclass) names1 <- paste("Class ",1:object$nclass)
		else names1 <- ''
		names2 <- names(object$patterns)[1:object$blocksize]
		names2 <- strsplit(names2,"\\.")
		x <- NULL
		for (i in 1:object$blocksize) {
			x <- c(x,names2[[i]][1])
		}
		if (object$blocksize==1) names2 <- ""
		else names2 <- x
		dimnames(out$lambdacoef) <- list(names1,names2)
		if (object$level2) {
			out$taucoef <- object$taucoef
			if (object$byclass) {
				names(out$taucoef) <- paste("Class ",1:object$nclass)
			} else {
				names(out$taucoef) <-''
			}
		}
	}
	class(out) <- "summary.randomLCA"
	out
}


print.summary.randomLCA <- function(x, ...)
{
	print(data.frame(Classes = x$nclass, AIC = x$AIC, BIC = x$BIC,
		log.Lik = c(x$log.Lik),row.names = " ") )
    cat("Class probabilities","\n")
	print(x$classp,digits=4)
	if (x$random)  cat("Conditional outcome probabilities","\n")
    else cat("Outcome probabilities","\n")
	print(x$outcomep,digits=4)
	if (x$random) {
		cat("Loadings","\n")
		if (length(x$lambdacoef)==1) cat(sprintf("%g\n",x$lambdacoef))
		else print(x$lambdacoef,digits=4)
		if (x$level2) {
			cat("Tau","\n")
			if (length(x$taucoef)==1) cat(sprintf("%g\n",x$taucoef))
			else print(x$taucoef,digits=4)
		}
	}
	invisible(x)
}