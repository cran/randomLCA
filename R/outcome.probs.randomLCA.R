# taken from confint
format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
	  "%")

outcome.probs <- function(object, ...)
	UseMethod("outcome.probs")

outcome.probs.randomLCA <-
function(object,level = 0.95,...) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
    if (object$random)
        stop("Not implemented for models with random effects yet.\n")
    out <- list()
# needs changes for random effects models
	nclass <- object$nclass
	outcomex <- object$fit$estimate[nclass:(length(object$fit$estimate))]
	outcomexse <- object$se[nclass:(length(object$fit$estimate))]
# transform using logistic to probabilities  
	if (object$probit) {
		outcomep <- pnorm(outcomex)
		outcomepl <- pnorm(outcomex+qnorm((1-level)/2)*outcomexse)
		outcomepu <- pnorm(outcomex+qnorm(1 - (1-level)/2)*outcomexse)
	} else {
		outcomep <- exp(outcomex)/(1+exp(outcomex))
		outcomepl <- exp(outcomex+qnorm((1-level)/2)*outcomexse)/(1+exp(outcomex+qnorm((1-level)/2)*outcomexse))
		outcomepu <- exp(outcomex+qnorm(1 - (1-level)/2)*outcomexse)/(1+exp(outcomex+qnorm(1 - (1-level)/2)*outcomexse))
	}
	
	outcomep <- t(matrix(outcomep,nrow=nclass))
	outcomepl <- t(matrix(outcomepl,nrow=nclass))
	outcomepu <- t(matrix(outcomepu,nrow=nclass))

	qnames <- format.perc(c((1-level)/2,1 - (1-level)/2),3)
	
	for (i in 1:nclass) {
		oneoutcome <- data.frame(outcomep[,i],outcomepl[,i],outcomepu[,i])
		row.names(oneoutcome) <- names(object$patterns)
		names(oneoutcome) <- c("Outcome p",qnames)
		out <- c(out,NULL)
		out[[i]] <- oneoutcome
	}
	class(out) <- "outcome.probs.randomLCA"
	out
}


print.outcome.probs.randomLCA <- function(x, ...)
{
	for (i in 1:length(x)) {
		cat("Class ",i,"\n")
		print(x[[i]])
	}
	invisible(x)
}

