`plot.randomLCA` <-
function(x,...,graphtype=ifelse(x$random,"marginal","conditional"),conditionalp=0.5,classhorizontal=TRUE) {
    if (!inherits(x, "randomLCA"))
        stop("Use only with 'randomLCA' xs.\n")
    if (missing(graphtype)) graphtype <- ifelse(x$random,"marginal","conditional")
    # calculate everything to be plotted
    graphdata <- NULL
    if (graphtype=="marginal") {
    	graphdata <- calc.marg.prob(x)
    	graphdata <- cbind(perc=rep(0,dim(graphdata)[1]),graphdata)
    }
    if (graphtype=="conditional")  graphdata <- calc.cond.prob(x,conditionalp)
    if (graphtype=="conditional2") graphdata <- calc.cond2.prob(x,conditionalp)
    
# set up the x axis labels    
#	thenames <- names(object$patterns)[1:object$blocksize]
#		names2 <- strsplit(names2,"\\.")
#		x <- NULL
#		for (i in 1:object$blocksize) {
#			x <- c(x,names2[[i]][1])
#		}
#		names2 <- x
#  ???? need to do something with blocksize
  if (x$level2) blocksize <- x$level2size
  else blocksize <- x$blocksize
# first decide if there are multiple blocks, otherwise plot all classes on one graph
    if (blocksize==dim(x$patterns)[2]) {
			if (graphtype=="marginal")
				print(xyplot(outcomep~outcome,group=class,data=graphdata,...))
			if ((graphtype=="conditional") || (graphtype=="conditional2")) {
				if (length(conditionalp)>1)
					print(xyplot(outcomep~outcome|perc,group=class,data=graphdata,...))
				else print(xyplot(outcomep~outcome,group=class,data=graphdata,...))
		}
    }
     else {
    	if (classhorizontal) {
			if (graphtype=="marginal")
				print(xyplot(outcomep~block|class,group=outcome,data=graphdata,...))
			if ((graphtype=="conditional") || (graphtype=="conditional2")) {
				if (length(conditionalp)>1)
					print(xyplot(outcomep~block|class*perc,group=outcome,data=graphdata,...))
				else print(xyplot(outcomep~block|class,group=outcome,data=graphdata,...))
			}
		}
		else  {
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

