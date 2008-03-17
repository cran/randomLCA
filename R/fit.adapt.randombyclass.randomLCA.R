`fit.adapt.randombyclass.randomLCA` <- function(patterns,freq,nclass,calcSE,initoutcomep,initclassp,initlambdacoef,gh,blocksize,probit,verbose) {

	nrepeats <- dim(patterns)[2]/blocksize
	nlevel1 <- dim(patterns)[2]
	nlevel2 <- length(freq)

	if (verbose) print("fit.random.randomLCA")

		# first rearrange the patterns
		npatterns <- 1-patterns
		# if an NA setting to zero in both patterns excludes from calculation
		myoutcomes <- ifelse(is.na(as.matrix(patterns)),0,as.matrix(patterns))
		mynoutcomes <- ifelse(is.na(as.matrix(npatterns)),0,as.matrix(npatterns))

		newmyoutcomes <- NULL
		for (i in 1:length(gh[,1])) newmyoutcomes <- rbind(newmyoutcomes,myoutcomes)
		newmynoutcomes <- NULL
		for (i in 1:length(gh[,1])) newmynoutcomes <- rbind(newmynoutcomes,mynoutcomes)

		calclikelihood <- function(classx,outcomex,lambdacoef,momentdata,gh) {
			#starttime <- proc.time()
			# turn classx into actual probabilities
			classp2 <- c(0,classx)       
			classp2 <- exp(classp2)/sum(exp(classp2))
			# calculate the likelihood over all classes and level 2 units
			ltotal <- NULL
			ll <- matrix(rep(NA,nclass*nlevel2),ncol=nclass)
			ll2 <- matrix(rep(NA,length(gh[,1])*nlevel2),nrow=nlevel2)
			level2p <- momentdata[,1]+momentdata[,2] %o% gh[,1]
			level2w <- log(sqrt(2*pi))+log(momentdata[,2])+
				matrix(rep(gh[,1]^2/2+log(gh[,2]),each=nlevel2),nrow=nlevel2)+
				matrix(dnorm(as.vector(level2p),log=TRUE),nrow=nlevel2)
			for (iclass in 1:nclass) {
				# outcome data
				myoutcomex <- outcomex[iclass,]
				# now start calculating for each level 2 quadrature point
				myoutcomex2 <- t(myoutcomex+t(as.vector(level2p) %o% rep(exp(lambdacoef[iclass,]),nrepeats)))
				if (probit) {
					lmyoutcomep <- pnorm(myoutcomex2,log.p=TRUE)
					nlmyoutcomep <- pnorm(-myoutcomex2,log.p=TRUE)
				}
				else {
					lmyoutcomep <- -log(1+exp(-myoutcomex2))
					nlmyoutcomep <- -log(1+exp(myoutcomex2))
				}
				ll2 <- matrix(rowSums(newmyoutcomes*lmyoutcomep+newmynoutcomes*nlmyoutcomep),ncol=length(gh[,1]))
				# final level 2 likelihoods
				ll[,iclass] <- log(rowSums(exp(ll2+level2w),na.rm=TRUE))
			}
			ill <- t(t(exp(ll))*classp2)
			ill2 <- rowSums(ill,na.rm=TRUE)
			ll <- sum(log(ill2)*freq,na.rm=TRUE)
			fitted <- ill2*sum(ifelse(apply(patterns,1,function(x) any(is.na(x))),0,freq))*
				ifelse(apply(patterns,1,function(x) any(is.na(x))),NA,1)
			classprob <- ill/ill2
			return(list(logl=ll,fitted=fitted,classprob=classprob))
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
				  return(ll)
				}
			  optim.fit <- nlm(loglik,x[length(x)],print.level=0,iterlim=1000,hessian=TRUE,outcomes=x[1:(length(x)-1)],gradtol=1.0e-7)
#  calculate se
			  return(c(beta=optim.fit$estimate[1],sebeta=sqrt(1/optim.fit$hessian)))
		}
		betas <- t(apply(cbind(patterns,momentdata[,1]),1,onerandom))
		return(betas)
	}
	
	adaptivefit <- function(classx,outcomex,lambdacoef,momentdata,calcSE,gh) {
	
	
		fitparams <- function(classx,outcomex,lambdacoef,
			momentdata,calcSE,gh,noiterations=10) {
			calcllfornlm <- function(params,momentdata,gh) {
				oneiteration <- calclikelihood(if (nclass==1) NULL else params[1:(nclass-1)],
					matrix(params[nclass:(nclass+nlevel1*nclass-1)],nrow=nclass),
					matrix(params[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*blocksize-1)],nrow=nclass),
					momentdata,gh)
				ll <- -oneiteration$logl
				ll
			}
			
			nlm1 <- nlm(calcllfornlm, c(classx, as.vector(outcomex), lambdacoef), iterlim = noiterations,
				print.level=ifelse(verbose,2,0),hessian=calcSE,
				check.analyticals = FALSE,momentdata=momentdata,gh=gh)
			return(list(logl=-nlm1$minimum,
				classx=(if (nclass==1) NULL else nlm1$estimate[1:(nclass-1)]),
				outcomex=matrix(nlm1$estimate[nclass:(nclass+nlevel1*nclass-1)],ncol=dim(patterns)[2]),
				lambdacoef=matrix(nlm1$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*blocksize-1)],nrow=nclass),
				nlm=nlm1))	
		}
	

		oneiteration <- calclikelihood(classx,outcomex,lambdacoef,
			momentdata,gh)
		currll <- oneiteration$logl
		if (verbose) cat('Initial ll',currll,"\n")
	# shift the quadrature points for the first time
		momentdata <- calcrandom(classx,outcomex,lambdacoef,momentdata)
		oneiteration <- calclikelihood(classx,outcomex,lambdacoef,
		momentdata,gh)
		currll <- oneiteration$logl
		if (verbose) cat("current ll",currll,"\n")		
		
		currll <- 1
		
		last2ll <- currll*2
		while(abs((last2ll-currll)/last2ll)>1.0e-6) {
			last2ll <- currll
			# need to do an optimisation on the other parameters
			fitresults <- fitparams(classx,outcomex,lambdacoef,
				momentdata,FALSE,gh)
			currll <- fitresults$logl
			outcomex <- fitresults$outcomex
			classx <- fitresults$classx
			lambdacoef <- fitresults$lambdacoef
			if (verbose) cat("current ll from optimisation",currll,"\n")		
			# shift the quadrature points again
			momentdata <- calcrandom(classx,outcomex,lambdacoef,momentdata)
			# fix up any strange se
			momentdata[,2] <- ifelse(is.finite(momentdata[,2]),momentdata[,2],1)
			oneiteration <- calclikelihood(classx,outcomex,lambdacoef,
			momentdata,gh)
			currll <- oneiteration$logl
			if (verbose) cat("current ll",currll,"\n")		
		}
		fitresults <- fitparams(classx,outcomex,lambdacoef,
				momentdata,calcSE,gh,noiterations=500)
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
			oneiteration <- calclikelihood(classx,outcomex,lambdacoef,
				momentdata,gh)
			currll <- oneiteration$logl
			if (verbose) cat('Initial ll',currll,"\n")
			lastll <- 2*currll
		# shift the quadrature points for the first time
			while (abs((lastll-currll)/lastll)>1.0e-6) {
				lastll <- currll
				momentdata <- calcrandom(classx,outcomex,lambdacoef,momentdata)
				oneiteration <- calclikelihood(classx,outcomex,lambdacoef,
				momentdata,gh)
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
 	
 	lambdacoef <- ifelse(lambdacoef < -3,-3,lambdacoef)
  	   
 	myfit <- adaptivefit(classx,outcomex,lambdacoef,momentdata,calcSE,gh)

	optim.fit <- myfit$nlm
	momentdata <- myfit$momentdata

	classx <- NULL
	if (nclass>1) classx <- optim.fit$estimate[1:(nclass-1)]
	outcomex <- matrix(optim.fit$estimate[nclass:(nclass+nclass*nlevel1-1)],ncol=dim(patterns)[2])
# transform using logistic to probabilities     
    lambdacoef <- matrix(optim.fit$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*blocksize-1)],nrow=nclass)
	
	final <- calclikelihood(classx,outcomex,lambdacoef,
			momentdata,gh)
			
	fitted <- final$fitted
	classprob <- final$classprob

# calculate the probabilities
# add extra column to classx
	classx <- c(0,classx)       
	 classp <- exp(classx)/apply(matrix(exp(classx),nrow=1),1,sum)

    if (probit) outcomep <- pnorm(outcomex)
    else outcomep <- exp(outcomex)/(1+exp(outcomex))
	
# extract the se
	if (!calcSE) separ <- rep(NA,length(optim.fit$estimate))
	else {
		s <- svd(optim.fit$hessian)
		separ <- sqrt(diag(s$v %*% diag(ifelse(s$d==0,NA,1/s$d)) %*% t(s$u)))
	}

# determine random effects

	ranef <- calcrandom(classx,outcomex,lambdacoef,momentdata)


	
    np <- length(optim.fit$estimate)
    aic <- 2*optim.fit$minimum+2*np
    bic <- 2*optim.fit$minimum+log(sum(freq))*np
    deviance <- 2*sum(ifelse(freq==0,0,freq*log(freq/fitted)))
    list(fit=optim.fit,nclass=nclass,classp=classp,outcomep=outcomep,lambdacoef=lambdacoef,se=separ,
    	aic=aic,bic=bic,np=np,log.Lik=-optim.fit$minimum,observed=freq,fitted=fitted,
    	deviance=deviance,ranef=ranef,classprob=classprob)
}
