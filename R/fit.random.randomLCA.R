`fit.random.randomLCA` <-
function(patterns,freq,nclass,calcSE,initoutcomep,initclassp,initlambdacoef,gh,blocksize,probit,verbose) {

	nrepeats <- dim(patterns)[2]/blocksize
	nlevel1 <- dim(patterns)[2]
	nlevel2 <- length(freq)

	noutcomes <- 1-patterns
	# if an NA setting to zero in both patterns excludes from calculation
	myoutcomes <- ifelse(is.na(as.matrix(patterns)),0,as.matrix(patterns))
	mynoutcomes <- ifelse(is.na(as.matrix(noutcomes)),0,as.matrix(noutcomes))
	
    loglik <- function(params) {
        if (nclass==1) classp <- 1
        else {
            classx <- params[1:(nclass-1)]
        	classx <- ifelse(abs(classx)>10,sign(classx)*10,classx)
# add extra column to classp
             classx <- c(0,classx)
# transform using logistic to probabilities		
       		classp <- exp(classx)/sum(exp(classx))
        }       
        outcomex <- matrix(params[nclass:(length(params)-blocksize)],ncol=nlevel1)
        outcomex <- ifelse(abs(outcomex)>10,sign(outcomex)*10,outcomex)
        lambdacoef <- params[(length(params)-blocksize+1):length(params)]
#		calculate probabilities for each random effect level
		ill <- matrix(rep(NA,nclass*length(freq)),ncol=nclass)
		probs <- matrix(rep(NA,length(gh[,1])*length(freq)),ncol=length(gh[,1]))
# calculate probabilities under each class
        for (i in 1:nclass) {
			newoutcomex <- outcomex[i,]+rep(lambdacoef,nrepeats) %o% gh[,1]
			if(probit) lnewoutcomep <- pnorm(newoutcomex,log.p=TRUE)
			else lnewoutcomep <- -log(1+exp(-newoutcomex))
			if(probit) nlnewoutcomep <- pnorm(-newoutcomex,log.p=TRUE)
			else nlnewoutcomep <- -log(1+exp(newoutcomex))
			probs <- exp(myoutcomes %*% lnewoutcomep+mynoutcomes %*% nlnewoutcomep)			
# multiply by class probabilities
			ill[,i] <- probs %*% as.vector(gh[,2])*classp[i]
        }
		ill2 <- rowSums(ill)
        ll <- -sum(log(ill2)*freq)
		attr(ll, "fitted") <- ill2*sum(ifelse(apply(patterns,1,function(x) any(is.na(x))),0,freq))*ifelse(apply(patterns,1,function(x) any(is.na(x))),NA,1)
		attr(ll, "classprob") <- ill/ill2
		return(ll)
    }
	
	calclikelihood <- function(classx,outcomex,lambdacoef) {
		-loglik(c(classx,as.vector(outcomex),lambdacoef))
	}

    if (nclass==1) classx <- NULL
    else  {
        classx <- rep(NA,nclass-1)
    	initclassp <- ifelse(initclassp<1.0e-4,1.0e-4,initclassp)       	
    	initclassp <- ifelse(initclassp>1-1.0e-4,1-1.0e-4,initclassp)
    	initclassp <- initclassp/sum(initclassp)
        for (i in 2:nclass) classx[i-1] <- log(initclassp[i]/initclassp[1])
    }

	initoutcomep <- ifelse(initoutcomep<1.0e-4,1.0e-4,initoutcomep)
	initoutcomep <- ifelse(initoutcomep>(1-1.0e-4),1-1.0e-4,initoutcomep)
    if (probit) outcomex <- qnorm(initoutcomep)
    else outcomex <- log(initoutcomep/(1-initoutcomep))

# choose among possible lambdacoef
 	if (missing(initlambdacoef) || is.null(initlambdacoef)) {
 		testlambdacoef <- 0
 		maxllambda <- NA
 		maxll <- -Inf
 		repeat {
			if (verbose) cat('trying lambdacoef ',testlambdacoef,"\n")
			lambdacoef <- rep(testlambdacoef,blocksize)
			currll <- calclikelihood(classx,outcomex,lambdacoef)
			if (verbose) cat("ll",currll,"\n")		
			# when the ll starts decreasing, give up
			if (currll < maxll) break()
			maxll <- currll
			maxllambda <- testlambdacoef
			testlambdacoef <- testlambdacoef+0.1
		}
		if (verbose) cat('using lambdacoef ',testlambdacoef,"\n")
		lambdacoef <- rep(maxllambda,blocksize)
 	}
 	else lambdacoef <- initlambdacoef

    optim.fit <- nlm(loglik,c(as.vector(classx),as.vector(outcomex),lambdacoef),hessian=calcSE,
    	print.level=ifelse(verbose,2,0),iterlim=1000)

    if (optim.fit$code >= 3)
    	warning("nlm exited with code ",optim.fit$code," .\n")
	if (nclass==1) classx <- NULL
	else classx <- optim.fit$estimate[1:(nclass-1)]
	outcomex <- matrix(optim.fit$estimate[nclass:(nclass+nlevel1*nclass-1)],ncol=nlevel1)
# transform using logistic to probabilities     
    lambdacoef <- optim.fit$estimate[(nclass+nlevel1*nclass):length(optim.fit$estimate)]
	
	final <- calclikelihood(classx,outcomex,lambdacoef)
			
	fitted <- attr(final,"fitted")
	classprob <- attr(final,"classprob")

# calculate the probabilities
    if (nclass==1) classp <- 1
    else {
# add extra column to classx
        classx <- c(0,classx)       
   		 classp <- exp(classx)/sum(exp(classx))
    }       
    if (probit) outcomep <- pnorm(outcomex)
    else outcomep <- exp(outcomex)/(1+exp(outcomex))
	
# extract the se
	if (!calcSE) separ <- rep(NA,length(optim.fit$estimate))
	else {
		s <- svd(optim.fit$hessian)
		separ <- sqrt(diag(s$v %*% diag(ifelse(s$d==0,NA,1/s$d)) %*% t(s$u)))
	}

# determine random effects
	calcrandom <- function() {
			
		onerandom <- function(x) {
	
				loglik <- function(beta) {
			# calculate probabilities under each class
						for (i in 1:nclass) {
			# calculate the outcome probabilities for this class and current random
							if (probit) outcomep <- pnorm(outcomex[i,]+rep(lambdacoef,nrepeats)*beta)
							else outcomep <- 1/(1+exp(-outcomex[i,]-rep(lambdacoef,nrepeats)*beta))
							oneprob <- t(apply(t(x)*outcomep+t(1-x)*(1-outcomep),2,prod,na.rm=TRUE))
			# multiply by class probabilities
							if (i==1) allprob <- oneprob*classp[i]
							else allprob <- allprob+oneprob*classp[i]
						}
					
					ll <- -(sum(log(allprob))+dnorm(beta,mean=0,sd=1,log=TRUE))
				  return(ll)
				}
			  beta <- 0
			  optim.fit <- nlm(loglik,beta,print.level=0,iterlim=1000,hessian=TRUE)
#  calculate se
			  return(c(beta=optim.fit$estimate[1],sebeta=sqrt(1/optim.fit$hessian)))
		}
		betas <- t(apply(patterns,1,onerandom))
		return(betas)
	}

	ranef <- calcrandom()


	
    np <- length(optim.fit$estimate)
    aic <- 2*optim.fit$minimum+2*np
    bic <- 2*optim.fit$minimum+log(sum(freq))*np
    deviance <- 2*sum(ifelse(freq==0,0,freq*log(freq/fitted)))
    list(fit=optim.fit,nclass=nclass,classp=classp,outcomep=outcomep,lambdacoef=lambdacoef,se=separ,
    	aic=aic,bic=bic,np=np,log.Lik=-optim.fit$minimum,observed=freq,fitted=fitted,
    	deviance=deviance,ranef=ranef,classprob=classprob)
}

