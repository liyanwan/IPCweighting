print.dIPCWSuperLearner <- function(x, ...) {
    cat("\nCall: ", deparse(x$call, width.cutoff = .9*getOption("width")), "\n\n", fill = getOption("width"))
    print(cbind(Risk = x$cvRisk, Coef = x$coef))
}

coef.dIPCWSuperLearner <- function(object, ...) {
    object$coef
}

