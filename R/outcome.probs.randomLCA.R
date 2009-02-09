# taken from confint
format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
	  "%")

outcome.probs <- function(object, ...)
	UseMethod("outcome.probs")

outcome.probs.randomLCA <-
function(object,level = 0.95, boot=FALSE, type="perc",R=ifelse(type=="perc",999,200),...) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
    if (object$random & !boot)
        stop("Not implemented for models with random effects. Use bootstrap option.\n")
    out <- list()
	nclass <- object$nclass
	if (boot) {
		dostatistic <- function(x,initmodel) {
	# fit the model using the simulated data
			onesim <- function () {
				sim <- randomLCA(x,nclass=initmodel$nclass,
					calcSE=FALSE,initmodel=initmodel,
					blocksize=initmodel$blocksize,notrials=1,random=initmodel$random,
					byclass=initmodel$byclass,quadpoints=initmodel$quadpoints,
					level2=initmodel$level2,probit=initmodel$probit,verbose=FALSE)
				if (sim$random) {
					marg.prob <- calc.marg.prob(sim)
					marg.prob <- as.vector(matrix(marg.prob$outcomep,
						nrow=nclass,byrow=TRUE))
# convert to appropriate scale
					if (sim$probit) estimate <- qnorm(marg.prob)
					else estimate <- log(marg.prob/(1-marg.prob))	
					} else  estimate <- sim$fit$estimate[nclass:(length(sim$fit$estimate))]
					return(estimate)
				}
			estimate <- tryCatch(onesim(),
				error=function(e) {
					warning(e)
					#browser()
					return(rep(NA,dim(initmodel$outcomep)[1]*dim(initmodel$outcomep)[2]))
				})					
			return(estimate)
		}
		gendata <- function(d, p) {
			myseed <- as.integer(runif(1, 0, .Machine$integer.max))
			simdata <- simulate(p, nsim=1,seed=myseed)[[1]]
			return(simdata)
		}
# expand data so one observation per row
		newdata <- object$patterns[rep(1:length(object$freq),object$freq),]
		newmodel <- randomLCA(newdata,nclass=object$nclass,calcSE=FALSE,initmodel=object,
				blocksize=object$blocksize,notrials=1,random=object$random,
				byclass=object$byclass,quadpoints=object$quadpoints,
				level2=object$level2,probit=object$probit)
		themle <- newmodel
		theboot <- boot(newdata,dostatistic,R=R,sim="parametric",
		ran.gen=gendata,
			mle=themle,
			initmodel=newmodel
			)
		failed <- sum(is.na(theboot$t[,1]))
		if (failed > 0) warning(sprintf("Convergence failed for %d models",failed))
		ci <- t(apply(as.matrix(1:length(theboot$t0)),1,function(x) {
				ci <- boot.ci(theboot, conf = level, type = type,index = x)
				if (type=="perc") return(c(ci$t0,ci$percent[4:5]))
				if (type=="normal") return(c(ci$t0,ci$normal[2:3]))
				if (type=="basic") return(c(ci$t0,ci$basic[4:5]))
		}))
        
		outcomex <- ci[,1]
		outcomexl <- ci[,2]
		outcomexu <- ci[,3]
		
# transform using logistic to probabilities  
		if (object$probit) {
			outcomep <- pnorm(outcomex)
			outcomepl <- pnorm(outcomexl)
			outcomepu <- pnorm(outcomexu)
		} else {
			outcomep <- exp(outcomex)/(1+exp(outcomex))
			outcomepl <- exp(outcomexl)/(1+exp(outcomexl))
			outcomepu <- exp(outcomexu)/(1+exp(outcomexu))
		}
	} else {
# needs changes for random effects models
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

