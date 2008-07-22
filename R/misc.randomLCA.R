logLik.randomLCA <- function(object, ...)
{
     val <- object$logLik
    attr(val, "df") <- object$np
    attr(val, "nobs") <- object$nobs
    class(val) <- "logLik"
    val
}

AIC.logLik <-
  ## BIC for logLik objects
  function(object, ...)
{
  -2 * (c(object) - attr(object, "df"))
}

BIC.logLik <-
  ## BIC for logLik objects
  function(object, ...)
{
  -2 * (c(object) - attr(object, "df") * log(attr(object, "nobs"))/2)
}


#AIC <-
#  ## Return the object's value of the Bayesian Information Criterion
#  function(object, ...) UseMethod("BIC")
#
#BIC <-
#  ## Return the object's value of the Bayesian Information Criterion
#  function(object, ...) UseMethod("BIC")

AIC.randomLCA <-
  ## AIC for various fitted objects
  function(object, ...)
{
	AIC(logLik(object))
}

BIC.randomLCA <-
  ## BIC for various fitted objects
  function(object, ...)
{
	BIC(logLik(object))
}
