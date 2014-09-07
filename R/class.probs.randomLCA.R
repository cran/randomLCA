class.probs <- function(object)
  UseMethod("class.probs")

class.probs.randomLCA <-
  function(object) {
    if (!inherits(object, "randomLCA"))
      stop("Use only with 'randomLCA' objects.\n")
  object$classp
}   