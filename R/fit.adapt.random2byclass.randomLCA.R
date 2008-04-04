fit.adapt.random2byclass.randomLCA <- function(outcomes,freq,nclass=2,initoutcomep,initclassp,initlambdacoef,initltaucoef,
    blocksize,calcSE=FALSE,gh,probit,verbose=FALSE) {

	if (probit) stop('not implemented')
	
    nlevel1 <- blocksize
    nlevel2 <- dim(outcomes)[2]/blocksize
    nlevel3 <- length(freq)

    # first rearrange the outcomes
      myoutcomes <- t(matrix(t(outcomes),nrow=nlevel1))
      mynoutcomes <- t(matrix(t(1-outcomes),nrow=nlevel1))
      
      myoutcomes <- ifelse(is.na(myoutcomes),0,myoutcomes)
      mynoutcomes <- ifelse(is.na(mynoutcomes),0,mynoutcomes)
      

    calclikelihood <- function(classx,outcomex,lambdacoef,ltaucoef,momentdata,gh,updatemoments=FALSE) {

		# turn classp into actual probabilities
		classp2 <- c(0,classx)       
		classp2 <- exp(classp2)/sum(exp(classp2))
		# rearrage moment data to allow extraction by class
		momentdata <- matrix(as.vector(momentdata),byrow=T,nrow=nclass)
		classdata <- cbind(momentdata,outcomex)
		# calculate the likelihood over all classes and level 3 units
        ltotal <- NULL
        ll <- matrix(rep(NA,nclass*nlevel3),ncol=nclass)
        ll2 <- matrix(rep(NA,length(gh[,1])*nlevel3*nlevel2),nrow=nlevel3*nlevel2)
        if (updatemoments) {
        # temporary for calculating moments by level 3 quadrature
            mye_2 <- matrix(rep(NA,length(gh[,1])*nlevel3*nlevel2),nrow=nlevel3*nlevel2)
            mye2_2 <- matrix(rep(NA,length(gh[,1])*nlevel3*nlevel2),nrow=nlevel3*nlevel2)
            mye23_2 <- matrix(rep(NA,length(gh[,1])*nlevel3*nlevel2),nrow=nlevel3*nlevel2)
        # level 2 moments after integration
            e_2 <- matrix(rep(NA,nclass*nlevel3*nlevel2),nrow=nlevel3*nlevel2)
            e2_2 <- matrix(rep(NA,nclass*nlevel3*nlevel2),nrow=nlevel3*nlevel2)
            e23_2 <- matrix(rep(NA,nclass*nlevel3*nlevel2),nrow=nlevel3*nlevel2)
        # level 3 moments after integration
            e <- matrix(rep(NA,nclass*nlevel3),nrow=nlevel3)
            e2 <- matrix(rep(NA,nclass*nlevel3),nrow=nlevel3)
        }

        level2ll <- matrix(rep(NA,length(gh[,1])*nlevel3*nlevel2),nrow=nlevel3*nlevel2)
        for (iclass in 1:nclass) {
            x <- classdata[iclass,]
            mymomentdata <- matrix(x[1:(nlevel3*(2+3*nlevel2))],nrow=nlevel3)
            level2moments <- t(matrix(t(mymomentdata[,3:(2+3*nlevel2)]),nrow=3))
            level3moments <- mymomentdata[,1:2]
            
            level3p <- level3moments[,1]+level3moments[,2] %o% gh[,1]
            level3w <- log(sqrt(2*pi))+log(level3moments[,2])+
                matrix(rep(gh[,1]^2/2+log(gh[,2]),each=nlevel3),nrow=nlevel3)+
                matrix(dnorm(as.vector(level3p),log=TRUE),nrow=nlevel3)
    # expand the level 3 quadrature points
            level3mu <- rep(mymomentdata[,1],each=nlevel2)
            rlevel3p <- matrix(rep(as.vector(level3p),each=nlevel2),nrow=nlevel2*nlevel3)
            rlevel3w <- matrix(rep(as.vector(level3w),each=nlevel2),nrow=nlevel2*nlevel3)
            # outcome data
            myoutcomex <- 
                t(matrix(rep(x[(nlevel3*(2+3*nlevel2)+1):(nlevel3*(2+3*nlevel2)+nlevel2*nlevel1)],nlevel3),
                nrow=nlevel1))
            # now start calculating for each level 3 quadrature point
            for (i3 in 1:length(gh[,1])) {
                # calculate location of level 2 quadrature points
                level2p <- level2moments[,1]+level2moments[,2] %o% gh[,1]
                level2p <- level2p-level2moments[,3]*(rlevel3p[,i3]-level3mu)
                level2w <- 0.5*(log(2)+log(pi))+log(level2moments[,2])+
                    matrix(rep(gh[,1]^2/2+log(gh[,2]),each=nlevel3*nlevel2),nrow=nlevel3*nlevel2)+
                    matrix(dnorm(as.vector(level2p),log=TRUE),nrow=nlevel3*nlevel2)
                # calculate for each level 2 quadrature point
                for (i2 in 1:length(gh[,1])) {
                    # calculate the outcome probabilities
                    myoutcomex2 <- myoutcomex+ (rlevel3p[,i3]+level2p[,i2]*exp(ltaucoef[iclass])) %o% lambdacoef[iclass,]
                    #myoutcomep <- 1/(1+exp(-myoutcomex2))
                    #level2ll[,i2] <- rowSums(log(myoutcomes*myoutcomep+mynoutcomes*(1-myoutcomep)))
					if(probit) lmyoutcomep <- pnorm(myoutcomex2,log.p=TRUE)
					else lmyoutcomep <- -log(1+exp(-myoutcomex2))
					if(probit) nlmyoutcomep <- pnorm(-myoutcomex2,log.p=TRUE)
					else nlmyoutcomep <- -log(1+exp(myoutcomex2))
                    level2ll[,i2] <- rowSums(myoutcomes*lmyoutcomep+mynoutcomes*nlmyoutcomep)
                }
                # calculate the total likelihood and moments for each level 2 unit
                ll2[,i3] <- log(rowSums(exp(level2w+level2ll)))
                if (updatemoments) {
                    mye_2[,i3] <- rowSums(level2p*exp(level2w+level2ll))/
                        exp(ll2[,i3])
                    mye2_2[,i3] <- rowSums(level2p^2*exp(level2w+level2ll))/
                        exp(ll2[,i3])
                    mye23_2[,i3] <- rowSums(rlevel3p[,i3]*level2p*exp(level2w+level2ll))/exp(ll2[,i3])
                }
            }
            # level 3 likelihoods by quadrature points
           # ll3 <- matrix(colSums(matrix(as.vector(ll2),nrow=nlevel2),na.rm=T),nrow=nlevel3)
            ll3 <- matrix(colSums(matrix(as.vector(ll2),nrow=nlevel2)),nrow=nlevel3)
            # final level 3 likelihoods
            ll[,iclass] <- log(rowSums(exp(ll3+level3w)))
            if (updatemoments) {
                # calculate level 3 moments
                e[,iclass] <- rowSums(level3p*exp(level3w+ll3))/
                    exp(ll[,iclass])
                e2[,iclass] <- rowSums((level3p-e[,iclass])^2*exp(level3w+ll3))/
                    exp(ll[,iclass])
                e2[,iclass] <- sqrt(e2[,iclass])
                # integrate level 2 results over level 3
                rll3 <- matrix(rep(as.vector(ll3),each=nlevel2),nrow=nlevel3*nlevel2)
                e_2[,iclass] <- rowSums(mye_2*exp(rll3+rlevel3w))/
                    exp(rep(ll[,iclass],each=nlevel2))
                e2_2[,iclass] <- rowSums(mye2_2*exp(rll3+rlevel3w))/
                    exp(rep(ll[,iclass],each=nlevel2))
                e2_2[,iclass] <- sqrt(e2_2[,iclass]-e_2[,iclass]^2)
                e23_2[,iclass] <- rowSums(mye23_2*exp(rll3+rlevel3w))/
                    exp(rep(ll[,iclass],each=nlevel2))
                e23_2[,iclass] <- -(e23_2[,iclass]-rep(e[,iclass],each=nlevel2)*e_2[,iclass])/
                    (rep(e2[,iclass]^2,each=nlevel2))
            }
        }
        # set up moments in correct form
        ltotal <- NULL
        if (updatemoments) {
            for (iclass in 1:nclass) {
                ltotal <- cbind(ltotal,e[,iclass],e2[,iclass])
                temp <- rbind(matrix(e_2[,iclass],ncol=nlevel2,byrow=T),
                                matrix(e2_2[,iclass],ncol=nlevel2,byrow=T),
                                matrix(e23_2[,iclass],ncol=nlevel2,byrow=T))
                temp <- matrix(temp,nrow=nlevel3)
                ltotal <- cbind(ltotal,temp)
            }
        }
        ill <- t(t(exp(ll))*classp2)
        ill2 <- log(rowSums(ill))
        totalll <- sum(ill2*freq)
        fitted <- exp(ill2)*sum(ifelse(apply(outcomes,1,function(x) any(is.na(x))),0,freq))*
            ifelse(apply(outcomes,1,function(x) any(is.na(x))),NA,1)
        classprob <- ill/exp(ill2)
        return(list(logl=totalll,moments=ltotal,fitted=fitted,classprob=classprob))
    }  # end of calclikelihood
    
    adaptivefit <- function(classx,outcomex,lambdacoef,ltaucoef,momentdata,gh) {
    
        fitparams <- function(classx,outcomex,lambdacoef,ltaucoef,
            momentdata,gh,calcSE,noiterations=10) {
                    
            calcllfornlm <- function(params,momentdata,gh) {  
                oneiteration <- calclikelihood(if (nclass==1) NULL else params[1:(nclass-1)],
                	matrix(params[nclass:(nclass+nlevel1*nlevel2*nclass-1)],nrow=nclass),
                    matrix(params[(nlevel1*nlevel2*nclass+nclass):(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1-1)],nrow=nclass),
                    params[(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1):(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1+nclass-1)],
                    momentdata,gh)
                return(-oneiteration$logl)
            }
                        
            nlm1 <- nlm(calcllfornlm, c(classx,as.vector(outcomex),lambdacoef,ltaucoef),
            	iterlim = noiterations,
                print.level=ifelse(verbose,2,0),
                check.analyticals = FALSE,hessian=calcSE,momentdata=momentdata,gh=gh)
            return(list(logl=-nlm1$minimum,
            	classx=if (nclass==1) NULL else nlm1$estimate[1:(nclass-1)],
            	outcomex=matrix(nlm1$estimate[nclass:(nclass+nlevel1*nlevel2*nclass-1)],nrow=nclass),
                lambdacoef=matrix(nlm1$estimate[(nlevel1*nlevel2*nclass+nclass):(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1-1)],nrow=nclass),
                ltaucoef=nlm1$estimate[(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1):(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1+nclass-1)],
                nlm=nlm1)) 
        }
    

       oneiteration <- calclikelihood(classx, outcomex, lambdacoef,ltaucoef,
            momentdata,gh,updatemoments=TRUE)
        currll <- oneiteration$logl
        if (verbose) cat('Initial ll',currll,"\n")
        lastll <- 2*currll
    # shift the quadrature points for the first time
        while (abs((lastll-currll)/lastll)>1.0e-6) {
            lastll <- currll
            momentdata <- oneiteration$moments
            oneiteration <- calclikelihood(classx,outcomex,lambdacoef,ltaucoef,
            momentdata,gh,updatemoments=TRUE)
            currll <- oneiteration$logl
            if (verbose) cat("current ll",currll,"\n")       
        }
                
		adaptive <- TRUE
		prevll <- -Inf
        while(adaptive) {
            # need to do an optimisation on the other parameters
            fitresults <- fitparams(classx,outcomex,lambdacoef,ltaucoef,
                momentdata,gh,calcSE=FALSE)
            currll <- fitresults$logl
            outcomex <- fitresults$outcomex
            classx <- fitresults$classx
            lambdacoef <- fitresults$lambdacoef
            ltaucoef <- fitresults$ltaucoef
            if (verbose) cat("current ll from optimisation",currll,"\n") 
            optll <- currll
            # shift the quadrature points again
            oneiteration <- calclikelihood(classx,outcomex,lambdacoef,ltaucoef,
                momentdata,gh,updatemoments=TRUE)
            currll <- oneiteration$logl
            lastll <- 2*currll
            while(abs((lastll-currll)/lastll)>1.0e-7) {
                lastll <- currll
                momentdata <- oneiteration$moments
                oneiteration <- calclikelihood(classx,outcomex,lambdacoef,ltaucoef,
                momentdata,gh,updatemoments=TRUE)
                currll <- oneiteration$logl
            if (verbose) cat("current ll",currll,"\n")       
            }
        	adaptive <- (abs((currll-optll)/currll)>1.0e-7) ||
        		(abs((currll-prevll)/currll)>1.0e-7)
        	if ((prevll-currll)/abs(currll) > 1.0e-4) stop("divergence - increase quadrature points")
        	prevll <- currll
       }
        fitresults <- fitparams(classx,outcomex,lambdacoef,ltaucoef,
                momentdata,gh,calcSE=calcSE,noiterations=500)
        return(list(nlm=fitresults$nlm,momentdata=momentdata))
    } # end adaptivefit

    # momentdata is level3 and level2
    # mu3,tau3,nlevel2*(mu2,taucoef,gamma2)
    # repeated for each class
    
        momentdata <- matrix(rep(c(rep(c(0,1),each=nlevel3),
                rep(rep(c(0,1,0),each=nlevel3),times=nlevel2)),times=nclass),
                nrow=nlevel3)
    
    if (nclass==1) classx <- NULL
    else  {
        classx <- rep(NA,nclass-1)
        initclassp <- ifelse(initclassp<1.0e-4,1.0e-4,initclassp)        
        initclassp <- ifelse(initclassp>(1.0-1.0e-4),1-1.0e-4,initclassp)          
        for (i in 2:nclass) classx[i-1] <- log(initclassp[i]/initclassp[1])
        }
        
    initoutcomep <- ifelse(initoutcomep<1.0e-4,1.0e-4,initoutcomep)
    initoutcomep <- ifelse(initoutcomep>(1.0-1.0e-4),1-1.0e-4,initoutcomep)
    if (probit) outcomex <- qnorm(initoutcomep)
    else outcomex <- log(initoutcomep/(1-initoutcomep))
    if (missing(initlambdacoef) || is.null(initlambdacoef)) lambdacoef <- matrix(rep(0,blocksize*nclass),nrow=nclass)
    else lambdacoef <- initlambdacoef

# choose among possible ltaucoef
 	if (missing(initltaucoef) || is.null(initltaucoef)) {
 		testltaucoef <- -3
 		maxltau <- NA
 		maxll <- -Inf
 		repeat {
			if (verbose) cat('trying ltaucoef ',testltaucoef,"\n")
			testltaucoef <- rep(testltaucoef,nclass)
			onelikelihood <- calclikelihood(classx,outcomex,lambdacoef,testltaucoef,momentdata,gh)
			currll <- onelikelihood$logl
			if (verbose) cat("ll",currll,"\n")		
			# when the ll starts decreasing, give up
			if (currll < maxll) break()
			maxll <- currll
			maxltau <- testltaucoef
			testltaucoef <- testltaucoef+0.1
		}
		if (verbose) cat('using ltaucoef ',testltaucoef,"\n")
		ltaucoef <- maxltau
 		ltaucoef <- rep(ltaucoef,nclass)
 	}
 	else ltaucoef <- initltaucoef
 	
        
    myfit <- adaptivefit(classx, outcomex, lambdacoef,ltaucoef,momentdata,gh)
    
    optim.fit <- myfit$nlm
    momentdata <- myfit$momentdata
    
# extract the se
    if (!calcSE) separ <- rep(NA,length(optim.fit$estimate))
    else {
		s <- svd(optim.fit$hessian)
		separ <- sqrt(diag(s$v %*% diag(ifelse(s$d==0,NA,1/s$d)) %*% t(s$u)))
    }
# calculate the probabilities
    if (nclass==1) classp <- 1
    else {
        classp <-optim.fit$estimate[1:(nclass-1)]
# add extra column to classp
        classp <- c(0,classp)       
    }       
    outcomep <- matrix(optim.fit$estimate[nclass:(nclass+nlevel1*nlevel2*nclass-1)],nrow=nclass)
# transform using logistic to probabilities     
    classp <- exp(classp)/sum(exp(classp))
    outcomep <- exp(outcomep)/(1+exp(outcomep))

	lambdacoef <- matrix(optim.fit$estimate[(nlevel1*nlevel2*nclass+nclass):(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1-1)],nrow=nclass)
    ltaucoef <- optim.fit$estimate[(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1):(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1+nclass-1)]

	calcrandom <- function() {
	
		outcomex <- log(outcomep/(1-outcomep))
			
		onerandom <- function(x) {
		
				loglik <- function(beta) {
			# calculate probabilities under each class
						for (i in 1:nclass) {
			# calculate the outcome probabilities for this class and current random
							outcomep <- 1/(1+exp(-outcomex[i,]-rep(lambdacoef[i,],nlevel2)*(beta[1]+exp(ltaucoef[i])*rep(beta[2:(1+nlevel2)],each=nlevel1))))
							oneprob <- t(apply(t(x)*outcomep+t(1-x)*(1-outcomep),2,prod,na.rm=TRUE))
			# multiply by class probabilities
							if (i==1) allprob <- oneprob*classp[i]
							else allprob <- allprob+oneprob*classp[i]
						}
					
					ll <- -(sum(log(allprob))+sum(dnorm(beta,mean=0,sd=1,log=TRUE)))
				  return(ll)
				}
			  beta <- rep(0,1+nlevel2)
			  optim.fit <- nlm(loglik,beta,print.level=0,iterlim=1000,hessian=TRUE,gradtol=1.0e-7)
				beta <- optim.fit$estimate
				sebeta <- sqrt(diag(solve(optim.fit$hessian)))
				checkx <- matrix(x,ncol=nlevel1,byrow=T)
				checkx <- apply(checkx,1,function(x) any(is.na(x)))
				checkx <- c(FALSE,checkx)
				beta[checkx] <- NA
				sebeta[checkx] <- NA
			  return(c(beta=beta,sebeta=sebeta))
		}
		betas <- t(apply(outcomes,1,onerandom))
		return(betas)
	}

    ranef <- calcrandom()

    final <- calclikelihood(if (nclass==1) NULL else optim.fit$estimate[1:(nclass-1)],
    				matrix(optim.fit$estimate[nclass:(nclass+nlevel1*nlevel2*nclass-1)],nrow=nclass),
                    matrix(optim.fit$estimate[(nlevel1*nlevel2*nclass+nclass):(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1-1)],nrow=nclass),
                    optim.fit$estimate[(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1):(nlevel1*nlevel2*nclass+nclass+nclass*nlevel1+nclass-1)],
                    momentdata,gh,updatemoments=FALSE)    
    fitted <- final$fitted
    classprob <- final$classprob
    
    np <- length(optim.fit$estimate)
    aic <- 2*optim.fit$minimum+2*np
    bic <- 2*optim.fit$minimum+log(sum(freq))*np
    deviance <- 2*sum(ifelse(freq==0,0,freq*log(freq/fitted)))
    list(fit=optim.fit,nclass=nclass,classp=classp,outcomep=outcomep,lambdacoef=lambdacoef,taucoef=exp(ltaucoef),se=separ,
        aic=aic,bic=bic,np=np,logl=-optim.fit$minimum,freq=freq,fitted=fitted,ranef=ranef
    ,classprob=classprob,deviance=deviance)
}
