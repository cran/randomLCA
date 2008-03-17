`fit.fixed.randomLCA` <-
function(patterns,freq,initoutcomep,initclassp,nclass,calcSE,verbose) {

# parameters
#   outcomes matrix of outcomes 0 or 1
#   freq vector of frequencies corresponding to each outcome combination
#   nclass number of classes
#   initoutcomep initial outcome probabilities
#   initclassp initial class probabilities
#   calcSE calculate standard errors ?
#   verbose print information about algorithm    

	nlevel1 <- dim(patterns)[2]
	
	npatterns <- 1-patterns
	# if an NA setting to zero in both patterns excludes from calculation
	mypatterns <- ifelse(is.na(as.matrix(patterns)),0,as.matrix(patterns))
	mynpatterns <- ifelse(is.na(as.matrix(npatterns)),0,as.matrix(npatterns))

    loglik <- function(params) {
        if (nclass==1) classp <- 1
        else {
            classx <- params[1:(nclass-1)]
# add extra column to classx
        	classx <- ifelse(abs(classx)>10,sign(classx)*10,classx)
             classx <- c(0,classx)
# transform using logistic to probabilities		
       		classp <- exp(classx)/apply(matrix(exp(classx),nrow=1),1,sum)
        }       
        outcomex <- matrix(params[nclass:length(params)],ncol=nlevel1)
        outcomex <- ifelse(abs(outcomex)>10,sign(outcomex)*10,outcomex)
		ill <- matrix(rep(NA,nclass*length(freq)),ncol=nclass)
# calculate probabilities under each class
        for (i in 1:nclass) {
# calculate the outcome probabilities for this class
			loutcomep <- -log(1+exp(-outcomex[i,]))
			nloutcomep <- -log(1+exp(outcomex[i,]))
			ill[,i] <- exp(mypatterns %*% loutcomep+ mynpatterns %*% nloutcomep)*classp[i]
# multiply by class probabilities
        }
		ill2 <- rowSums(ill)
        ll <- -sum(log(ill2)*freq)
		attr(ll, "fitted") <- ill2*sum(ifelse(apply(patterns,1,function(x) any(is.na(x))),0,freq))*ifelse(apply(patterns,1,function(x) any(is.na(x))),NA,1)
		attr(ll, "classprob") <- ill/ill2		
	  return(ll)
    }

	if (missing(initclassp))  initclassp <- runif(nclass)
	initclassp <- ifelse(initclassp<1.0e-3,1.0e-3,initclassp)       	
	initclassp <- ifelse(initclassp>1-1.0e-3,1-1.0e-3,initclassp)       	
	classp <- initclassp/sum(initclassp)
	
    if (missing(initoutcomep)) initoutcomep <- runif(nclass*nlevel1)
	initoutcomep <- ifelse(initoutcomep<1.0e-5,1.0e-5,initoutcomep)
	initoutcomep <- ifelse(initoutcomep>1-1.0e-5,1-1.0e-5,initoutcomep)
	outcomep <- matrix(initoutcomep,nrow=nclass)

# now do the em algorithm

	ill <- matrix(rep(NA,nclass*length(freq)),ncol=nclass)
	oldll <- NA
	emit <- 0
	repeat {
		# calculate probabilities under each class
		for (i in 1:nclass) {
			# calculate the outcome probabilities for this class
			ill[,i]  <- t(apply(t(patterns)*outcomep[i,]+t(1-patterns)*(1-outcomep[i,]),2,
							prod,na.rm=TRUE))*classp[i]
		}
		ill2 <- apply(ill,1,sum,na.rm=TRUE)
		ll <- sum(log(ill2)*freq,na.rm=TRUE)

		if(is.na(oldll)) oldll <- 2*ll
		if (abs((ll-oldll)/ll) < 1.0e-9) break()
		oldll <- ll
		
		# estimated posterior probability
		classprob <- ill/ill2
		
		# new classp
		classp <- apply(classprob,2,weighted.mean,w=freq)
		
		#new outcome probabilities
		for (i in 1:nclass) {
			outcomep[i,] <- apply(patterns*classprob[,i],2,weighted.mean,w=freq,na.rm=TRUE)/
								apply(ifelse(is.na(patterns),NA,1)*classprob[,i],2,weighted.mean,w=freq,na.rm=TRUE)
		}
		emit <- emit+1
		if (((emit %% 100)==0) & verbose) cat('iteration ',emit,' logl ',ll,'\n')
	}
	fitted <- ill2*sum(ifelse(apply(patterns,1,function(x) any(is.na(x))),0,freq))*ifelse(apply(patterns,1,function(x) any(is.na(x))),NA,1)
	# if SE required then use quasi-Newton
	if (calcSE) {
		outcomex <- log(outcomep/(1-outcomep))
		outcomex <- ifelse(abs(outcomex) > 6,6*sign(outcomex),outcomex)
		
		if (nclass==1) classx <- NULL
		else  {
			classx <- rep(NA,nclass-1)
			for (i in 2:nclass) classx[i-1] <- log(classp[i]/classp[1])
			classx <- ifelse(abs(classx) > 6,6*sign(classx),classx)
		}
	 
		  optim.fit <- nlm(loglik,c(as.vector(classx),as.vector(outcomex)),hessian=calcSE,print.level=ifelse(verbose,2,0),
			gradtol=1.0e-6,iterlim=1000)
		if (optim.fit$code >= 3)
			warning("nlm exited with code ",optim.fit$code," .\n")
		s <- svd(optim.fit$hessian)
		separ <- sqrt(diag(s$v %*% diag(ifelse(s$d==0,NA,1/s$d)) %*% t(s$u)))
# calculate the probabilities
		if (nclass==1) classx <- 0
		else {
			classx <- optim.fit$estimate[1:(nclass-1)]
# add extra column to classx
			classx <- c(0,classx)       
		}       
		outcomep <- matrix(optim.fit$estimate[nclass:(length(optim.fit$estimate))],ncol=nlevel1)
# transform using logistic to probabilities     
		classp <- exp(classx)/apply(matrix(exp(classx),nrow=1),1,sum)
		outcomep <- exp(outcomep)/(1+exp(outcomep))
	
		final <- loglik(optim.fit$estimate)
		fitted <- attr(final,"fitted")
		classprob <- attr(final,"classprob")
		ll <- -optim.fit$minimum
   }
    else {
    	optim.fit <- NULL
    	separ <- rep(NA,length(c(as.vector(classp),as.vector(outcomep))))
    }
    

    if (verbose) {
		print("results")
		print(classp)
		print(outcomep)
    }
    np <- (nclass-1)+nclass*nlevel1
    aic <- -2*ll+2*np
    bic <- -2*ll+log(sum(freq))*np
    deviance <- 2*sum(ifelse(freq==0,0,freq*log(freq/fitted)))
    list(fit=optim.fit,nclass=nclass,classp=classp,outcomep=outcomep,se=separ,
    	aic=aic,bic=bic,np=np,log.Lik=ll,observed=freq,fitted=fitted,
    	deviance=deviance,classprob=classprob)
}

