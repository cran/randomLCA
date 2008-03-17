`summary.randomLCA` <-
function(object,...) {
    if (!inherits(object, "randomLCA"))
        stop("Use only with 'randomLCA' objects.\n")
    print("LogL")
    print(object$log.Lik)
    print("Class probabilities")
    print(object$classp)
    print("Outcome Probabilities")
    print(object$outcomep)
    print("Loadings")
    print(object$lambdacoef)
    out <- list()
 out$log.Lik <- object$log.Lik
 out$AIC <- object$aic
 out$BIC <- object$bic
 class(out) <- "summ.randomLCA"
 invisible(out)
}

