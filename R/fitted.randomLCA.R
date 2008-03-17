`fitted.randomLCA` <-
function(object,...) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
     return(object$fitted)
 }

