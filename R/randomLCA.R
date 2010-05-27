`randomLCA` <-
function(patterns,freq,nclass=2,calcSE=TRUE,initmodel=NULL,blocksize=1,notrials=20,
	random=FALSE,byclass=FALSE,quadpoints=21,level2=FALSE,probit=FALSE,
	verbose=FALSE,seed = as.integer(runif(1, 0, .Machine$integer.max))) {
    set.seed(seed)
    if (quadpoints > 190)
        stop("Maximum of 190 quadrature points\n")
	cl <- match.call()
    # check that patterns are either 0 or 1
	if (any(apply(as.matrix(patterns),2,function(x) any((x!=0)&(x!=1)&!is.na(x)))))
		stop("patterns must consist of either 0 or 1")    
	# check that patterns doesn't contain column which is all missing
	if (any(apply(as.matrix(patterns),2,function(x) all(is.na(x)))))
		stop("patterns cannot contain columns consisting entirely of missing")
	# check that freq are all >= 0
	if (!missing(freq)) {
		if (any(freq<0))
			stop("frequencies must be greater than or equal to zero")
	}
	if (random & ((dim(patterns)[2] %% blocksize)!=0))
		stop("outcomes must be a multiple of blocksize")
	# if no frequencies given, then assume that the data needs to be summaried
	if (missing(freq)) {
		pats <- apply(patterns, 1, function(x) {paste(ifelse(is.na(x),"N",x),collapse="")})
		tpats <- table(pats)
		freq <- as.numeric(tpats)
		newpatterns <- unlist(strsplit(names(tpats),split=""))
		newpatterns <- ifelse(newpatterns=="N",NA_character_,newpatterns)
		newpatterns <- as.data.frame(matrix(as.numeric(newpatterns),byrow=TRUE,ncol=dim(patterns)[2]))
		if (is.null(names(patterns))) names(newpatterns) <- paste("X",1:dim(patterns)[2],sep="")
		else names(newpatterns) <- names(patterns)
		patterns <- newpatterns
	}
	else {
	# check that freq doesn't contain missing
		if (any(is.na(freq))) stop("freq cannot contain missing values")
	}
	if (missing(initmodel)) {
		initmodel <- bestlca(patterns,freq=freq,nclass=nclass,
		calcSE=(calcSE & !random),notrials=notrials,probit=probit,verbose=verbose)
		initmodel$nclass <- nclass
		initmodel$random <- FALSE
		initmodel$level2 <- FALSE
		initmodel$byclass <- FALSE
		initmodel$blocksize <- 1
	}
	else {
		if (initmodel$nclass != nclass)
			stop("Initial model must have same number of classes.\n")
		if (dim(initmodel$patterns)[2] != dim(patterns)[2])
			stop("Initial model must have same number of outcomes.\n")
		if (initmodel$random & !random)
			stop("Initial model random and model to be fitted not.\n")
		if (initmodel$random) {
			if ((initmodel$blocksize!=1) & (initmodel$blocksize!=blocksize))
				stop("Initial model must have same number or 1 for blocksize.\n")
			if (initmodel$byclass & !byclass)
				stop("Initial model by class and model to be fitted not.\n")
		}
		if (!random) initmodel <- fit.fixed.randomLCA(patterns,freq=freq,initoutcomep=initmodel$outcomep,initclassp=initmodel$classp,nclass=nclass,calcSE=calcSE,probit=probit,verbose=verbose)
	}
	if (!random) fit <- initmodel
	else {
		# sort out the initial lambdacoef
		if (!is.null(initmodel$lambdacoef)) {
			initlambdacoef <- as.vector(initmodel$lambdacoef)
			if (initmodel$blocksize != blocksize) initlambdacoef <- rep(initlambdacoef,each=blocksize)
			if (byclass) {
				if (!initmodel$byclass) {
					initlambdacoef <- rep(initlambdacoef,each=nclass)
				}
				initlambdacoef <- matrix(initlambdacoef,ncol=blocksize)
			}			
		}
		else {
			initlambdacoef <- NULL
		}
		if (level2) {
			if (!is.null(initmodel$taucoef)) {
				initltaucoef <- log(initmodel$taucoef)
				if (byclass) {
					if (initmodel$byclass != byclass) initltaucoef <- rep(initltaucoef,each=nclass)
				}
			}
			else initltaucoef <- NULL
		}
		if (!byclass) {
			if (level2) fit <- 	fit.adapt.random2.randomLCA(patterns,freq=freq,
					nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
					initclassp=initmodel$classp,initlambdacoef=initlambdacoef,
					initltaucoef=initltaucoef,
					gh=norm.gauss.hermite(quadpoints),blocksize=blocksize,
					probit=probit,verbose=verbose)
   			else {
				fit <- fit.adapt.random.randomLCA(patterns,freq=freq,nclass=nclass,
					calcSE=calcSE,initoutcomep=initmodel$outcomep,
					initclassp=initmodel$classp,initlambdacoef=initlambdacoef,
					gh=norm.gauss.hermite(quadpoints),
					blocksize=blocksize,probit=probit,verbose=verbose)
				}
			}
		else {
			if (level2) fit <- 	fit.adapt.random2byclass.randomLCA(patterns,freq=freq,
					nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
					initclassp=initmodel$classp,initlambdacoef=initlambdacoef,
					initltaucoef=initltaucoef,
					gh=norm.gauss.hermite(quadpoints),blocksize=blocksize,
					probit=probit,verbose=verbose)
   			else {
				fit <- fit.adapt.randombyclass.randomLCA(patterns,freq=freq,
					nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
					initclassp=initmodel$classp,initlambdacoef=initlambdacoef,
					gh=norm.gauss.hermite(quadpoints),
					blocksize=blocksize,probit=probit,verbose=verbose)
			}
		}
	}
	fit$call <- cl
	fit$nclass <- nclass
	fit$random <- random
	fit$level2 <- level2
	fit$byclass <- byclass
	fit$probit <- probit
	fit$quadpoints <- quadpoints
	fit$blocksize <- blocksize
	fit$patterns <- patterns
	fit$notrials <- notrials
	fit$freq <- freq
	class(fit) <- "randomLCA"
	return(fit)
}

