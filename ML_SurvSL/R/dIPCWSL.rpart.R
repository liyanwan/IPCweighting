#' Wrapper function for recursive partitioning and regression trees
#'
#' @param Y Outcome variable
#' @param X Covariate dataframe
#' @param newX Dataframe to predict the outcome
#' @param time_point cut-off time
#' @param family "binomial" for binary classification
#' @param id Optional id to group observations from the same unit (not used currently)
#' @param num_intervals The number of equal time intervals into which the time_point should be divided
#' @param cp Complexity parameter
#' @param minsplit Minimum number of observations required to attempt a split
#' @param ... Additional ignored arguments.
#' #' @seealso \code{\link{SL.rpart}} \code{\link[glmnet]{predict.dIPCWSL.rpart}}
#'   \code{\link[glmnet]{rpart}}

dIPCWSL.rpart <- function(Y, X, newX, time_point, family, id,
num_intervals = 1, cp = 0.01, minsplit = 20, xval = 0L, maxdepth = 30, minbucket = round(minsplit/3),...) {
  .SL.require('rpart')
  interval = seq(0, time_point, length.out = num_intervals + 1)
  expanded_train_data_with_weight = bin_combined_ipcw_separate(interval = interval, X = X, Y = Y,
                                                      time_point = time_point,  return_type = "dataset")
  varNames = colnames(X)
  valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
  filtered_train_E = expanded_train_data_with_weight$E[valid_indices]
  filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
  formula = as.formula(paste("E ~", paste(varNames, collapse = " + ")))
  model_data = expanded_train_data_with_weight[valid_indices,]
    fit.rpart <- rpart::rpart(formula, data = model_data, 
                              control = rpart::rpart.control(cp = cp, minsplit = minsplit, xval = xval, maxdepth = maxdepth, minbucket = minbucket), 
                              method = "class", weights = filtered_weights)
    pred_onebin <- predict(fit.rpart, newdata = newX)[, 2]
    pred <- 1 - (1 - pred_onebin)^num_intervals
    fit <- list(object = fit.rpart, num_intervals = num_intervals)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("dIPCWSL.rpart")
    return(out)
}

#' Prediction function for an dIPCWSL.glmnet object
#'
#' @param object Result object from dIPCWSL.rpart
#' @param newdata Dataframe or matrix that will generate predictions
#' @param family "binomial" for binary classification
#' @param ... Any additional arguments (not used).

predict.dIPCWSL.rpart <- function(object, newdata, family, ...) {
  .SL.require('rpart')
  if(family$family=="binomial") {
    pred_onebin <- predict(object$object, newdata = newdata)[, 2]
    num_intervals = object$num_intervals
    pred <-1 - (1 - pred_onebin)^num_intervals
  }
    return(pred)
}

