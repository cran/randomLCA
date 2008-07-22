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

	patterns <- as.matrix(patterns)
	mode(patterns) <- "double"
	
	nlevel1 <- dim(patterns)[2]

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
		ill <- matrix(rep(NA,nclass*length(freq)),ncol=nclass)
# calculate probabilities under each class
        for (i in 1:nclass) {
# calculate the outcome probabilities for this class
			outcomep <- 1/(1+exp(-outcomex[i,]))
			ill[,i] <- .Call("bernoulliprob",patterns,outcomep)*classp[i]
# multiply by class probabilities
        }
		ill2 <- rowSums(ill)
        ll <- -sum(log(ill2)*freq)
 	  return(ll)
    }
	
    calcfitted <- function(params) {
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
		ill <- matrix(rep(NA,nclass*length(freq)),ncol=nclass)
# calculate probabilities under each class
        for (i in 1:nclass) {
# calculate the outcome probabilities for this class
			outcomep <- 1/(1+exp(-outcomex[i,]))
			ill[,i] <- .Call("bernoulliprob",patterns,outcomep)*classp[i]
# multiply by class probabilities
        }
		ill2 <- rowSums(ill)
        ll <- -sum(log(ill2)*freq)
		fitted <- ill2*sum(ifelse(apply(patterns,1,function(x) any(is.na(x))),0,freq))*ifelse(apply(patterns,1,function(x) any(is.na(x))),NA,1)
		classprob <- ill/ill2
	  	return(list(fitted=fitted,classprob=classprob))
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
	x <- .Call("lcemalgorithm",patterns,outcomep,classp,as.numeric(freq),verbose)	
	ll <- x[[1]]
	outcomep <- x[[2]]
	classp <- x[[3]]

	outcomex <- log(outcomep/(1-outcomep))
	
	if (nclass==1) classx <- NULL
	else  {
		classx <- rep(NA,nclass-1)
		for (i in 2:nclass) classx[i-1] <- log(classp[i]/classp[1])
	}
 
	  optim.fit <- nlm(loglik,c(as.vector(classx),as.vector(outcomex)),hessian=calcSE,print.level=ifelse(verbose,2,0),
		gradtol=1.0e-6,iterlim=1000)
	if (optim.fit$code >= 3)
		warning("nlm exited with code ",optim.fit$code," .\n")
	if (calcSE) {
		s <- svd(optim.fit$hessian)
		separ <- sqrt(diag(s$v %*% diag(ifelse(s$d==0,NA,1/s$d)) %*% t(s$u)))
	} else separ <- rep(NA,length(c(as.vector(classp),as.vector(outcomep))));
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

	final <- calcfitted(optim.fit$estimate)
	fitted <- final$fitted
	classprob <- final$classprob
	ll <- -optim.fit$minimum
	if (verbose) {
		print("results")
		print(classp)
		print(outcomep)
    }
    np <- (nclass-1)+nclass*nlevel1
    nobs <- sum(freq)
    deviance <- 2*sum(ifelse(freq==0,0,freq*log(freq/fitted)))
    list(fit=optim.fit,nclass=nclass,classp=classp,outcomep=outcomep,se=separ,
    	np=np,nobs=nobs,logLik=ll,observed=freq,fitted=fitted,
    	deviance=deviance,classprob=classprob)
}

