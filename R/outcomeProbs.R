# taken from confint
format_perc <- function(probs, digits)
  ## Not yet exported, maybe useful in other contexts:
  ## quantile.default() sometimes uses a version of it
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
        "%")

outcomeProbs <- function(object,level = 0.95, boot=FALSE, type="perc",R=999,scale=c("prob","raw"),
                         cores = max(detectCores() %/% 2, 1),...)
  UseMethod("outcomeProbs")

outcomeProbs.randomLCA <-
  function(object,level = 0.95, boot=FALSE, type="perc",R=999,scale=c("prob","raw"),
           cores = max(detectCores() %/% 2, 1),...) {
    if (!inherits(object, "randomLCA"))
      stop("Use only with 'randomLCA' objects.\n")
    if (object$random & !boot)
      stop("Not implemented for models with random effects. Use bootstrap option.\n")
    if (object$random & boot & !(type %in% c("perc","norm")))
      stop("Only bootstrap of type perc or normal available")
    out <- list()
    nclass <- object$nclass
    if (boot) {
      dostatistic <- function(x,initmodel) {
        # fit the model using the simulated data
        onesim <- function () {
          #browser()
          sim <- refit(initmodel,newpatterns=x,useinit=TRUE)
          if (is.null(sim)) {
            estimate <- rep(NA,dim(initmodel$outcomep)[1]*dim(initmodel$outcomep)[2])
          } else {
          if (sim$random) {
            #print(sim$outcomep)
            marg.prob <- calcMargProb(sim)
            marg.prob <- as.vector(matrix(marg.prob$outcomep,
                                          nrow=nclass,byrow=TRUE))
            # convert to appropriate scale
            if (sim$probit) estimate <- qnorm(marg.prob)
            else estimate <- log(marg.prob/(1-marg.prob))	
          } else  estimate <- sim$fit$estimate[nclass:(length(sim$fit$estimate))]
          }
          return(estimate)
        }
        estimate <- tryCatch(onesim(),
                             error=function(e) {
                               if (cores==1) warning(e)
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
      newmodel <- refit(object,newpatterns=newdata,useinit=TRUE)
      themle <- newmodel
       if (cores>1) {
        if(.Platform$OS.type=="unix") parallel <- "multicore"
        else parallel <- "snow"
      } else parallel <- "no"
      theboot <- boot(newdata,dostatistic,R=R,sim="parametric",
                      ran.gen=gendata,
                      mle=themle,
                      initmodel=newmodel,
                      parallel=parallel,
                      ncpus = cores
      )
      failed <- sum(is.na(theboot$t[,1]))
      if (failed > 0) warning(sprintf("Convergence failed for %d models",failed))
       ci <- t(apply(as.matrix(1:length(theboot$t0)),1,function(x) {
        ci <- boot.ci(theboot, conf = level, type = type,index = x)
        if (type=="perc") return(c(ci$t0,ci$percent[4:5]))
        if (type=="norm") return(c(ci$t0,ci$normal[2:3]))
      }))
      outcomex <- ci[,1]
      outcomexl <- ci[,2]
      outcomexu <- ci[,3]
      
      if (scale[1]=="raw") {
        outcomep <- outcomex
        outcomepl <- outcomexl
        outcomepu <- outcomexu
      } else {
        # transform using logistic to probabilities
        if (object$probit) {
          outcomep <- pnorm(outcomex)
          outcomepl <- pnorm(outcomexl)
          outcomepu <- pnorm(outcomexu)
        } else {
          outcomep <- 1.0/(1+exp(-outcomex))
          outcomepl <- 1.0/(1+exp(-outcomexl))
          outcomepu <- 1.0/(1+exp(-outcomexu))
        }
      }
    } else {
      outcomex <- object$fit$estimate[nclass:(length(object$fit$estimate))]
      outcomexse <- object$se[nclass:(length(object$fit$estimate))]
      # transform using logistic or probit to probabilities  
      if (scale[1]=="raw") {
        outcomep <- outcomex
        outcomepl <- outcomex+qnorm((1-level)/2)*outcomexse
        outcomepu <- outcomex+qnorm(1 - (1-level)/2)*outcomexse
      } else {
        if (object$probit) {
          outcomep <- pnorm(outcomex)
          outcomepl <- pnorm(outcomex+qnorm((1-level)/2)*outcomexse)
          outcomepu <- pnorm(outcomex+qnorm(1 - (1-level)/2)*outcomexse)
        } else {
          outcomep <- 1.0/(1.0+exp(-outcomex))
          outcomepl <- 1.0/(1.0+exp(-outcomex-qnorm((1-level)/2)*outcomexse))
          outcomepu <- 1.0/(1.0+exp(-outcomex-qnorm(1 - (1-level)/2)*outcomexse))
        }
      }
    }
    
    outcomep <- t(matrix(outcomep,nrow=nclass))
    outcomepl <- t(matrix(outcomepl,nrow=nclass))
    outcomepu <- t(matrix(outcomepu,nrow=nclass))
    
    qnames <- format_perc(c((1-level)/2,1 - (1-level)/2),3)
    
    for (i in 1:nclass) {
      oneoutcome <- data.frame(outcomep[,i],outcomepl[,i],outcomepu[,i])
      row.names(oneoutcome) <- names(object$patterns)
      names(oneoutcome) <- c("Outcome p",qnames)
      out <- c(out,NULL)
      out[[i]] <- oneoutcome
    }
    class(out) <- "outcomeProbs.randomLCA"
    out
  }


print.outcomeProbs.randomLCA <- function(x, ...)
{
  for (i in 1:length(x)) {
    cat("Class ",i,"\n")
    print(x[[i]])
  }
  invisible(x)
}

