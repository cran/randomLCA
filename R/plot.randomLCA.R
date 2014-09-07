`plot.randomLCA` <-
function(x,...,graphtype=c("conditional","marginal","conditional2","both"),conditionalp=0.5,classhorizontal=TRUE) {
    if (!inherits(x, "randomLCA"))
        stop("Use only with 'randomLCA' xs.\n")
    if (missing(graphtype)) graphtype <- "conditional"
    # ???? not certain if this is correct way
    if (length(graphtype) > 1) graphtype <- graphtype[1]
    if ((graphtype=="both") & (length(conditionalp) > 1))
        stop("Only single conditional probability allowed for both graphs.\n") 
    # calculate everything to be plotted
    graphdata <- NULL
    if ((graphtype=="marginal") | (graphtype=="both")) {
    	graphdata <- calc.marg.prob(x)
    	graphdata <- cbind(perc=rep(0,dim(graphdata)[1]),graphdata)
    	graphdata$graphtype <- "Marginal"
    }
    if ((graphtype=="conditional") | (graphtype=="both")) {
    	conddata <- calc.cond.prob(x,conditionalp)
    	conddata$graphtype <- "Conditional"
    	graphdata <- rbind(graphdata,conddata)
    }
    if (graphtype=="conditional2") {
    	conddata <- calc.cond2.prob(x,conditionalp)
    	conddata$graphtype <- "Conditional2"
    	graphdata <- rbind(graphdata,conddata)
    }
    
# set up the x axis labels    
#	thenames <- names(object$patterns)[1:object$blocksize]
#		names2 <- strsplit(names2,"\\.")
#		x <- NULL
#		for (i in 1:object$blocksize) {
#			x <- c(x,names2[[i]][1])
#		}
#		names2 <- x
#  ???? need to do something with blocksize
	graphdata$graphtype <- factor(graphdata$graphtype)
  if (x$level2) blocksize <- x$level2size
  else blocksize <- x$blocksize
# first decide if there are multiple blocks, otherwise plot all classes on one graph
    if ((blocksize==1) | (blocksize==dim(x$patterns)[2])) {
    	if (blocksize==1) {
			if (graphtype=="both")
				print(xyplot(outcomep~block|graphtype,group=class,data=graphdata,...))
			if (graphtype=="marginal")
				print(xyplot(outcomep~block,group=class,data=graphdata,...))
			if (graphtype=="conditional") {
				if (length(conditionalp)>1) print(xyplot(outcomep~block|perc,group=class,data=graphdata,...))
				else print(xyplot(outcomep~block,group=class,data=graphdata,...))
			}
		}
		else {
			if (graphtype=="both")
				print(xyplot(outcomep~outcome|graphtype,group=class,data=graphdata,...))
			if (graphtype=="marginal")
				print(xyplot(outcomep~outcome,group=class,data=graphdata,...))
			if ((graphtype=="conditional") || (graphtype=="conditional2")) {
				if (length(conditionalp)>1)
					print(xyplot(outcomep~outcome|perc,group=class,data=graphdata,...))
				else print(xyplot(outcomep~outcome,group=class,data=graphdata,...))
			}
		}

     }
     else {
    	if (classhorizontal) {
			if (graphtype=="both")
				print(xyplot(outcomep~block|class*graphtype,group=outcome,data=graphdata,...))
			if (graphtype=="marginal")
				print(xyplot(outcomep~block|class,group=outcome,data=graphdata,...))
			if ((graphtype=="conditional") || (graphtype=="conditional2")) {
				if (length(conditionalp)>1)
					print(xyplot(outcomep~block|class*perc,group=outcome,data=graphdata,...))
				else print(xyplot(outcomep~block|class,group=outcome,data=graphdata,...))
			}
		}
		else  {
			if (graphtype=="both")
				print(xyplot(outcomep~block|graphtype*class,group=outcome,data=graphdata,...))
			if (graphtype=="marginal")
				print(xyplot(outcomep~block|class,group=outcome,data=graphdata,...))
			if ((graphtype=="conditional") || (graphtype=="conditional2")) {
				if (length(conditionalp)>1)
					print(xyplot(outcomep~block|perc*class,group=outcome,data=graphdata,...))
				else print(xyplot(outcomep~block|class,group=outcome,data=graphdata,...))
			}
		}
    }
}

