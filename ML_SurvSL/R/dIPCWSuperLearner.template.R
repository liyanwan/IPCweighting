# 
dIPCWSuperLearner.template <- function(Y, X, newX, id, ...) {
		# pred returns predicted responses (on the scale of the outcome)
		pred <- numeric()
		# fit <- list(object = )
		fit <- vector("list", length=0)
		class(fit) <- c("dIPCWSL.template")
		out <- list(pred = pred, fit = fit)
		return(out)
}

# 
predict.dIPCWSuperLearner.template <- function(object, newdata, X = NULL, Y = NULL,...) {
	# insert prediction function
	pred <- numeric()
	return(pred)
}

write.dIPCWSuperLearner.template <- function(file = '', ...) {
  cat('dIPCWSuperLearner.template <- function(Y, X, newX, id, ...) {\n  # load required packages\n  # require(\'pkg\')\n \n 		# insert estimation and prediction function \n # pred is the predicted responses for newX (on the scale of the outcome)\n  pred <- numeric()\n  # fit returns all objects needed for predict.dIPCWSuperLearner.template\n  fit <- list(object = )\n  # declare class of fit for predict.dIPCWSuperLearner.template\n  class(fit) <- \'dIPCWSuperLearner.template\'\n  # return a list with pred and fit\n  out <- list(pred = pred, fit = fit)\n  return(out)\n}', file = file, ...)
}