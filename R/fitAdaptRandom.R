`fitAdaptRandom` <- function(patterns,freq,nclass,calcSE,justEM,
                             initoutcomep,initclassp,initlambdacoef,initmomentdata,gh,constload,blocksize,probit,byclass,qniterations,
                             penalty,EMtol,verbose) {
  
  # parameters
  #   outcomes matrix of outcomes 0 or 1
  #   freq vector of frequencies corresponding to each outcome combination
  #   nclass number of classes
  #   initoutcomep initial outcome probabilities
  #   initclassp initial class probabilities
  #   initlambdacoef initial lambda coefficient
  #   constload are the loadings constant or different for each outcome
  #   calcSE calculate standard errors ?
  #   gh matrix of gauss-hermite coefficients first column positions, second columns weights
  #   probit use probit transform rather than logitic to obtain outcome probabilities
  #   verbose print information about algorithm    
  
  
  calcrandom <- function(classx,outcomex,lambdacoef,momentdata) {
    
    classx <- c(0,classx)
    classp <- exp(classx)/sum(exp(classx))
    
    onerandom <- function(x) {
      
      loglik <- function(beta,outcomes) {
        # calculate probabilities under each class
        for (i in 1:nclass) {
          # calculate the outcome probabilities for this class and current random
          if (byclass) {
            if (probit) outcomep <- pnorm(outcomex[i,]+rep(lambdacoef[i,],nrepeats)*beta)
            else outcomep <- 1/(1+exp(-outcomex[i,]-rep(lambdacoef[i,],nrepeats)*beta))
            if (probit) noutcomep <- pnorm(-outcomex[i,]-rep(lambdacoef[i,],nrepeats)*beta)
            else noutcomep <- 1/(1+exp(outcomex[i,]+rep(lambdacoef[i,],nrepeats)*beta))
          } else {
            if (probit) outcomep <- pnorm(outcomex[i,]+rep(lambdacoef,nrepeats)*beta)
            else outcomep <- 1/(1+exp(-outcomex[i,]-rep(lambdacoef,nrepeats)*beta))
            if (probit) noutcomep <- pnorm(-outcomex[i,]-rep(lambdacoef,nrepeats)*beta)
            else noutcomep <- 1/(1+exp(outcomex[i,]+rep(lambdacoef,nrepeats)*beta))
          }
          oneprob <- exp(sum(ifelse(outcomes==1,log(outcomep),0)+ifelse(outcomes==0,log(noutcomep),0),na.rm=TRUE))
          # multiply by class probabilities
          if (i==1) allprob <- oneprob*classp[i]
          else allprob <- allprob+oneprob*classp[i]
        }
        ll <- -(log(allprob)+dnorm(beta,mean=0,sd=1,log=TRUE))
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
  
  
  calclikelihood <- function(classx,outcomex,lambdacoef,momentdata,gh,patterns,freq,calcfitted=FALSE,zprop=NULL) {
    #starttime <- proc.time()
    # turn classx into actual probabilities
    classp2 <- c(0,classx)
    classp2 <- exp(classp2)/sum(exp(classp2))
    ill <- matrix(rep(NA,nclass*length(freq)),ncol=nclass)
    for (iclass in 1:nclass) {
      if (byclass) {
        #browser()
        ill[,iclass] <- .Call("bernoulliprobrandom",patterns,outcomex[iclass,],lambdacoef[iclass,],
                              gh,momentdata,probit)				} else {
                                ill[,iclass] <- .Call("bernoulliprobrandom",patterns,outcomex[iclass,],lambdacoef,
                                                      gh,momentdata,probit)
                              }
    }
    # if zprop not supplied then we have the usual maximum likelihood
    if (is.null(zprop)) {
      for (iclass in 1:nclass) ill[,iclass] <- ill[,iclass]+log(classp2[iclass])
      maxll <- rowMaxs(ill, value=TRUE)
      ll <- sum((maxll+log(rowSums(exp(ill-maxll))))*freq)
      ill2 <- rowSums(exp(ill))
    } else {
      # otherwise calculate the complete data maximum likelihood for the em algorithm
      maxll <- rowMaxs(ill+log(zprop), value=TRUE)
      ll <- sum((maxll+log(rowSums(exp(ill-maxll))))*freq)
      ill2 <- rowSums(exp(ill))
      ll <- sum(log(ill2)*freq,na.rm=TRUE)
    }
    # penalise extreme outcome probabilities
    if (probit) {
      outcomep <- pnorm(as.vector(outcomex))
      noutcomep <- pnorm(as.vector(-outcomex))
    } else {
      outcomep <- as.vector(1/(1+exp(-outcomex)))
      noutcomep <- as.vector(1/(1+exp(outcomex)))
    }
    #			pen <- dbeta(outcomep,1+penalty,1+penalty,log=TRUE)
    #			print(c(ll,sum(pen)))
    penll <- ll+penalty/(nclass*2)*sum(log(outcomep))+penalty/(nclass*2)*sum(log(noutcomep))
    #print(c(ll,penll,outcomep,noutcomep))
    if (is.nan(ll) || is.infinite(ll)) ll <- -1.0*.Machine$double.xmax
    if (calcfitted) {
      # do this again in case we are using likelihood for em
      ill2 <- rowSums(exp(ill),na.rm=TRUE)
      fitted <- ill2*sum(ifelse(apply(patterns,1,function(x) any(is.na(x))),0,freq))*
        ifelse(apply(patterns,1,function(x) any(is.na(x))),NA,1)
      classprob <- exp(ill)/ill2
      return(list(logLik=ll,penlogLik=penll,fitted=fitted,classprob=classprob))
    } else return(list(logLik=ll,penlogLik=penll))
  }  # end of calclikelihood
  
  fitparams <- function(classx,outcomex,lambdacoef,
                        momentdata,calcSE,gh,patterns,noiterations=qniterations,zprop=NULL) {
    calcllfornlm <- function(params,momentdata,gh,patterns,zprop) {
      oneiteration <- calclikelihood(if (nclass==1) NULL else params[1:(nclass-1)],
                                     matrix(params[nclass:(nclass+nlevel1*nclass-1)],nrow=nclass),
                                     if (byclass) matrix(params[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*lambdasize-1)],nrow=nclass)
                                     else params[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+lambdasize-1)],
                                     momentdata,gh,patterns,freq,zprop=zprop)
      ll <- -oneiteration$penlogLik
      if (is.nan(ll) || is.infinite(ll)) ll <- .Machine$double.xmax
      ll
    }
    
    nlm1 <- nlm(calcllfornlm, c(classx, as.vector(outcomex), lambdacoef), iterlim = noiterations,
                print.level=ifelse(verbose,2,0),hessian=calcSE,stepmax=1,
                check.analyticals = FALSE,momentdata=momentdata,gh=gh,patterns=patterns,zprop=zprop)
    return(list(penlogLik=-nlm1$minimum,
                classx=(if (nclass==1) NULL else nlm1$estimate[1:(nclass-1)]),
                outcomex=matrix(nlm1$estimate[nclass:(nclass+nlevel1*nclass-1)],nrow=nclass),
                lambdacoef=(if (byclass) matrix(nlm1$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+nclass*lambdasize-1)],nrow=nclass)
                            else lambdacoef=nlm1$estimate[(nlevel1*nclass+nclass):(nlevel1*nclass+nclass+lambdasize-1)]),
                nlm=nlm1))
  }
  
  #perform model fitting
  # momentdata is level2
  # mu2,sigma2,
  
  tryCatch({
    
    if (verbose) print("fit.random.randomLCA")
    
    patterns <- as.matrix(patterns)
    mode(patterns) <- "integer"
    
    nrepeats <- ifelse(constload,dim(patterns)[2],1)
    lambdasize <- ifelse(constload,1,min(dim(patterns)[2],blocksize))
    nlevel1 <- dim(patterns)[2]
    nlevel2 <- length(freq)
    
    if (nclass==1) classx <- NULL
    else {
      if (is.null(initclassp)) {
        classx <- rnorm(nclass-1)
        classx <- ifelse(classx > 2,2,classx)
        classx <- ifelse(classx < -2,-2,classx)
      } else  {
        classx <- rep(NA,nclass-1)
        initclassp <- ifelse(initclassp < 1.0e-3, 1.0e-3, initclassp)
        initclassp <-
          ifelse(initclassp > 1 - 1.0e-3, 1 - 1.0e-3, initclassp)
        initclassp <- initclassp/sum(initclassp)
        for (i in 2:nclass) classx[i-1] <- log(initclassp[i]/initclassp[1])
      }
    }
    
    if (is.null(initoutcomep)) {
      outcomex <- matrix(rnorm(nclass*dim(patterns)[2]),nrow=nclass)
      outcomex <- ifelse(outcomex > 2,2,outcomex)
      outcomex <- ifelse(outcomex < -2,-2,outcomex)
    }
    else {
      initoutcomep <- ifelse(initoutcomep<1.0e-4,1.0e-4,initoutcomep)
      initoutcomep <- ifelse(initoutcomep>(1-1.0e-4),1-1.0e-4,initoutcomep)
      if (probit) outcomex <- qnorm(initoutcomep)
      else outcomex <- log(initoutcomep/(1-initoutcomep))
    }
    
    if (is.null(initlambdacoef)) {
      if (constload) nlambda <- 1
      else nlambda <- min(dim(patterns)[2],blocksize)
      initlambdacoef <- 2*rbinom(nlambda,1,0.5)-1
      if (mean(initlambdacoef) < 0.0) initlambdacoef <- -initlambdacoef
      if (byclass) initlambdacoef <- matrix(rep(initlambdacoef,nclass),nrow=nclass,byrow=TRUE)
    }
    lambdacoef <- initlambdacoef
    
    if (is.null(initmomentdata)) momentdata <- matrix(rep(c(0,1),each=nlevel2),nrow=nlevel2)
    else momentdata <- initmomentdata
    
    momentdata <-  calcrandom(classx,outcomex,lambdacoef,momentdata)
    initlikelihood <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns,freq,calcfitted=TRUE,zprop=NULL)
    
    lastll <- -Inf
    currll <- initlikelihood$penlogLik
    
    prop <- initlikelihood$zprop
    
    
    #browser()
    
    while (abs((lastll-currll)/currll) > EMtol) {
      
      fitted.model <- fitparams(classx,outcomex,lambdacoef,
                                momentdata,calcSE+FALSE,gh,patterns,noiterations=qniterations,zprop=initlikelihood$zprop)
      
      
      classx <- fitted.model$classx
      outcomex <- fitted.model$outcomex
      lambdacoef <- fitted.model$lambdacoef
      momentdata <-  calcrandom(classx,outcomex,lambdacoef,momentdata)
      
      currlikelihood <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns,freq,calcfitted=TRUE,zprop=NULL)
      
      prop <- currlikelihood$zprop
      
      lastll <- currll
      currll <- currlikelihood$penlogLik 
      
      #print(c(currll,lastll))
    }
    
    #browser()
    
    if (!justEM) {
      
      lastll <- -Inf
      
      while (abs((lastll-currll)/currll) > 1.0e-7) {
        fitted.model <- fitparams(classx,outcomex,lambdacoef,
                                  momentdata,calcSE,gh,patterns,noiterations=qniterations,zprop=NULL)
        
        classx <- fitted.model$classx
        outcomex <- fitted.model$outcomex
        lambdacoef <- fitted.model$lambdacoef
        momentdata <-  calcrandom(classx,outcomex,lambdacoef,momentdata)

        currlikelihood <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns,freq,calcfitted=TRUE,zprop=NULL)
        
        lastll <- currll
        currll <- currlikelihood$penlogLik 
      }
      
      fitted.model <- fitparams(classx,outcomex,lambdacoef,
                                momentdata,calcSE,gh,patterns,noiterations=1000,zprop=NULL)
      
      classx <- fitted.model$classx
      outcomex <- fitted.model$outcomex
      lambdacoef <- fitted.model$lambdacoef
      momentdata <-  calcrandom(classx,outcomex,lambdacoef,momentdata)
      
      currlikelihood <- calclikelihood(classx,outcomex,lambdacoef,momentdata,gh,patterns,freq,calcfitted=TRUE,zprop=NULL)
      
    }
    
    
    # calculate the probabilities
    # add extra column to classx
    classx <- c(0,classx)
    classp <- exp(classx)/apply(matrix(exp(classx),nrow=1),1,sum)
    
    if (probit) outcomep <- pnorm(outcomex)
    else outcomep <- exp(outcomex)/(1+exp(outcomex))
    
    # extract the se
    if (!calcSE) separ <- rep(NA,length(fitted.model$nlm$estimate))
    else {
      s <- svd(fitted.model$nlm$hessian)
      separ <- diag(s$v %*% diag(1/s$d) %*% t(s$u))
      separ[!is.nan(separ) & (separ>=0.0)] <- sqrt(separ[!is.nan(separ) & (separ>=0.0)])
      separ[is.nan(separ) | (separ<0.0)] <- NA
    }
    
    # determine random effects
    
    ranef <- calcrandom(classx,outcomex,lambdacoef,momentdata)
    
    np <- length(fitted.model$nlm$estimate)
    nobs <- sum(freq)
    deviance <- 2*sum(ifelse(freq==0,0,freq*log(freq/currlikelihood$fitted)))
    #browser()
    return(list(fit=fitted.model$nlm,nclass=nclass,classp=classp,outcomep=outcomep,lambdacoef=lambdacoef,
                se=separ,momentdata=momentdata,
                np=np,nobs=nobs,logLik=currlikelihood$logLik,penlogLik=currlikelihood$penlogLik,observed=freq,
                fitted= currlikelihood$fitted,
                deviance=deviance,ranef=ranef,classprob=currlikelihood$classprob))
  },
  error=function(e) return(NULL)
  )
  
  
}
