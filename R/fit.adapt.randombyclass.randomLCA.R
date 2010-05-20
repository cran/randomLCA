`fit.adapt.randombyclass.randomLCA` <- function(patterns,freq,nclass,calcSE,initoutcomep,initclassp,initlambdacoef,gh,blocksize,probit,verbose) {

	patterns <- as.matrix(patterns)
	mode(patterns) <- "double"

	nrepeats <- dim(patterns)[2]/blocksize
	nlevel1 <- dim(patterns)[2]
	nlevel2 <- length(freq)

	if (verbose) print("fit.random.randomLCA")

		calclikelihood <- function(classx,outcomex,lambdacoef,momentdata,gh,patterns,calcfitted=FALSE) {
			#starttime <- proc.time()
			# turn classx into actual probabilities
			classp2 <- c(0,classx)       
			classp2 <- exp(classp2)/sum(exp(classp2))
			ill <- matrix(rep(NA,nclass*length(freq)),ncol=nclass)
			for (iclass in 1:nclass) {
				ill[,iclass] <- .Call("bernoulliprobrandom",patterns,outcomex[iclass,],lambdacoef[iclass,],
					gh,momentdata,probit)*classp2[iclass]
				#browser()
			}
			ill2 <- rowSums(ill,na.rm=TRUE)
			ll <- sum(log(ill2)*freq,na.rm=TRUE)
			if (is.nan(ll) || is.infinite(ll)) ll <- -1.0*.Machine$double.xmax
			#print(ll)
			if (calcfitted) {
				fitted <- ill2*sum(ifelse(apply(patterns,1,function(x) any(is.na(x))),0,freq))*
					ifelse(apply(patterns,1,function(x) any(is.na(x))),NA,1)
				classprob <- ill/ill2
				return(list(logl=ll,fitted=fitted,classprob=classprob))
			} else return(list(logl=ll))
		}  # end of calclikelihood

	calcrandom <- function(classx,outcomex,lambdacoef,momentdata) {

		classx <- c(0,classx)       
		classp <- exp(classx)/sum(exp(classx))
			
		onerandom <- function(x) {
	
				loglik <- function(beta,outcomes) {
			# calculate probabilities under each class
						for (i in 1:nclass) {
			# calculate the outcome probabilities for this class and current random
							if (probit) outcomep <- pnorm(outcomex[i,]+rep(lambdacoef[i,],nrepeats)*beta)
							else outcomep <- 1/(1+exp(-outcomex[i,]-rep(lambdacoef[i,],nrepeats)*beta))
							oneprob <- exp(sum(outcomes*log(outcomep)+(1-outcomes)*log(1-outcomep),na.rm=TRUE))
			# multiply by class probabilities
							if (i==1) allprob <- oneprob*classp[i]
							else allprob <- allprob+oneprob*classp[i]
						}
					
					ll <- -(sum(log(allprob))+dnorm(beta,mean=0,sd=1,log=TRUE))
					if (is.nan(ll) || is.infinite(ll)) ll <- .Machine$double.xmax
				  return(ll)
				}
			  optim.fit <- nlm(loglik,x[length(x)],print.level=0,iterlim=1000,hessian=TRUE,outcomes=x[1:(length(x)-1)],gradtol=1.0e-7)
#  calculate se
			  return(c(beta=optim.fit$estimate[1],sebeta=sqrt(1/optim.fit$hessian)))
		}
		betas <- t(apply(cbind(patterns,momentdata[,1]),1,onerandom))
		return(betas)
	}

	adaptivefit <- function(classx,outcomex,lambdacoef,calcSE,momentdata,gh,patterns) {
	
	
		fitparams <- function(classx,outcomex,lambdacoef,
			momentdata,calcSE,gh,patterns,noiterations=10) {
			calcllfornlm <- function(params,momentdata,gh,patterns) {
				oneiteration <- calclikelihood(if (nclass==1) NULL else params[1:(nclass-1)],
					matrix(params[nclass:(nclass+nlevel1*nclass-1)],nrow=nclass),
					matrix(params[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*blocksize-1)],nrow=nclass),
					momentdata,gh,patterns)
				ll <- -oneiteration$logl
				if (is.nan(ll) || is.infinite(ll)) ll <- .Machine$double.xmax
				ll
			}
			
			nlm1 <- nlm(calcllfornlm, c(classx, as.vector(outcomex), lambdacoef), iterlim = noiterations,
				print.level=ifelse(verbose,2,0),hessian=calcSE,
				check.analyticals = FALSE,momentdata=momentdata,gh=gh,patterns=patterns)
			return(list(logl=-nlm1$minimum,
				classx=(if (nclass==1) NULL else nlm1$estimate[1:(nclass-1)]),
				outcomex=matrix(nlm1$estimate[nclass:(nclass+nlevel1*nclass-1)],ncol=dim(patterns)[2]),
				lambdacoef=matrix(nlm1$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*blocksize-1)],nrow=nclass),
				nlm=nlm1))	
		}
	

		oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns)
		currll <- oneiteration$logl
		if (verbose) cat('Initial ll',currll,"\n")
	# shift the quadrature points for the first time
		momentdata <- calcrandom(classx,outcomex,lambdacoef,momentdata)
		oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns)
		currll <- oneiteration$logl
		if (verbose) cat("current ll",currll,"\n")		
		
		adaptive <- TRUE
		prevll <- -Inf
 		nadaptive <- 0
       while(adaptive) {
			# need to do an optimisation on the other parameters
			fitresults <- fitparams(classx,outcomex,lambdacoef,momentdata,FALSE,gh,patterns)
			currll <- fitresults$logl
			outcomex <- fitresults$outcomex
			classx <- fitresults$classx
			lambdacoef <- fitresults$lambdacoef
			if (verbose) cat("current ll from optimisation",currll,"\n")		
			# shift the quadrature points again
			momentdata <- calcrandom(classx,outcomex,lambdacoef,momentdata)
			oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns)
			# check if moving quadrature points has changed likelihood
        	adaptive <- (abs((oneiteration$logl-currll)/oneiteration$logl)>1.0e-7) ||
        		(abs((oneiteration$logl-prevll)/oneiteration$logl)>1.0e-7)
			currll <- oneiteration$logl
			if (verbose) cat("current ll",currll,"\n")
        	if ((prevll-currll)/abs(currll) > 1.0e-4) stop("divergence - increase quadrature points")
        	nadaptive <- nadaptive+1
        	if (nadaptive > 200) stop("too many adaptive iterations - increase quadrature points")
        	prevll <- currll
		}
		fitresults <- fitparams(classx,outcomex,lambdacoef,momentdata,calcSE,gh,patterns,noiterations=500)
		return(list(nlm=fitresults$nlm,momentdata=momentdata))
	} # end adaptivefit

	# momentdata is level2
	# mu2,lambda2,
	
	
    if (nclass==1) classx <- NULL
    else  {
        classx <- rep(NA,nclass-1)
    	initclassp <- ifelse(initclassp==0.0,1.0e-10,initclassp)       	
    	initclassp <- ifelse(initclassp==1.0,1-1.0e-10,initclassp)       	
        for (i in 2:nclass) classx[i-1] <- log(initclassp[i]/initclassp[1])
    }

	initoutcomep <- ifelse(initoutcomep<1.0e-3,1.0e-3,initoutcomep)
	initoutcomep <- ifelse(initoutcomep>(1-1.0e-3),1-1.0e-3,initoutcomep)
    if (probit) outcomex <- qnorm(initoutcomep)
    else outcomex <- log(initoutcomep/(1-initoutcomep))

	momentdata <- matrix(rep(c(0,1),each=nlevel2),nrow=nlevel2)

# choose among possible lambdacoef
 	if (missing(initlambdacoef) || is.null(initlambdacoef)) {
		lastmomentdata <- momentdata
 		testlambdacoef <- 0
 		maxllambda <- NA
 		maxll <- -Inf
 		repeat {
			if (verbose) cat('trying lambdacoef ',testlambdacoef,"\n")
			lambdacoef <- matrix(rep(testlambdacoef,nclass*blocksize),nrow=nclass)
			oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns)
			currll <- oneiteration$logl
			if (verbose) cat('Initial ll',currll,"\n")
			lastll <- 2*currll
		# shift the quadrature points for the first time
			while (abs((lastll-currll)/lastll)>1.0e-6) {
				lastll <- currll
				momentdata <- calcrandom(classx,outcomex,lambdacoef,momentdata)
				oneiteration <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns)
				currll <- oneiteration$logl
				if (verbose) cat("current ll",currll,"\n")		
			}
			momentdata <- calcrandom(classx,outcomex,lambdacoef,momentdata)
			# when the ll starts decreasing, give up
			if (currll < maxll) break()
			maxll <- currll
			maxllambda <- testlambdacoef
			lastmomentdata <- momentdata
			testlambdacoef <- testlambdacoef+0.1
		}
		if (verbose) cat('using lambdacoef ',maxllambda,"\n")
		lambdacoef <- matrix(rep(maxllambda,nclass*blocksize),nrow=nclass)
		momentdata <- lastmomentdata
 	}
 	else lambdacoef <- initlambdacoef
  	   
 	myfit <- adaptivefit(classx,outcomex,lambdacoef,calcSE,momentdata,gh,patterns)

	optim.fit <- myfit$nlm
	momentdata <- myfit$momentdata

	classx <- NULL
	if (nclass>1) classx <- optim.fit$estimate[1:(nclass-1)]
	outcomex <- matrix(optim.fit$estimate[nclass:(nclass+nclass*nlevel1-1)],ncol=dim(patterns)[2])
# transform using logistic to probabilities     
    lambdacoef <- matrix(optim.fit$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*blocksize-1)],nrow=nclass)
	
	final <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns,calcfitted=TRUE)
			
	fitted <- final$fitted
	classprob <- final$classprob

# calculate the probabilities
# add extra column to classx
	classx <- c(0,classx)       
	 classp <- exp(classx)/apply(matrix(exp(classx),nrow=1),1,sum)

    if (probit) outcomep <- pnorm(outcomex)
    else outcomep <- 1/(1+exp(-outcomex))
	
# extract the se
	if (!calcSE) separ <- rep(NA,length(optim.fit$estimate))
	else {
		s <- svd(optim.fit$hessian)
		separ <- sqrt(diag(s$v %*% diag(1/s$d) %*% t(s$u)))
		separ[!is.finite(separ)] <- NA
	}

# determine random effects

	ranef <- calcrandom(classx,outcomex,lambdacoef,momentdata)


	
    np <- length(optim.fit$estimate)
    nobs <- sum(freq)
    deviance <- 2*sum(ifelse(freq==0,0,freq*log(freq/fitted)))
    list(fit=optim.fit,nclass=nclass,classp=classp,outcomep=outcomep,lambdacoef=lambdacoef,
    	se=separ,
    	np=np,nobs=nobs,logLik=-optim.fit$minimum,observed=freq,fitted=fitted,
    	deviance=deviance,ranef=ranef,classprob=classprob)
}
