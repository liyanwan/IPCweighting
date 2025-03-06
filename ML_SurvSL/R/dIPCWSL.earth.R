#' Wrapper function using multivariate Adaptive Regression Splines
#'
#' @param Y Outcome variable
#' @param X Covariate dataframe
#' @param newX Dataframe to predict the outcome
#' @param time_point cut-off time
#' @param family "binomial" for binary classification.
#' @param id Optional id to group observations from the same unit (not used
#'   currently).
#' @param num_intervals The number of equal time intervals into which the time_point should be divided.
#' @param ... Additional ignored arguments.
#' 
dIPCWSL.earth <- function(Y, X, newX, time_point, family, id,
num_intervals = 1, degree = 2, nprune = 4,
penalty = if(degree>1) 3 else 2, nk = max(21, 2*ncol(X) + 1), pmethod = "backward", nfold = 0, ncross = 1, minspan = 0, endspan = 0,...) {
  .SL.require('earth')
  interval <- seq(0, time_point, length.out = num_intervals + 1)
  expanded_train_data_with_weight <- bin_combined_ipcw_separate(interval = interval, X = X, Y = Y,
                                                       time_point = time_point,  return_type = "dataset")
  varNames = colnames(X)
  train_weight_X <- expanded_train_data_with_weight[varNames]
  valid_indices <- which(expanded_train_data_with_weight$IPCW > 0)
  filtered_train_X <- train_weight_X[valid_indices, ]
  filtered_train_E <- expanded_train_data_with_weight$E[valid_indices]
  filtered_weights <- expanded_train_data_with_weight$IPCW[valid_indices]
  if(family$family == "binomial") {
    fit.earth <- earth::earth(x = filtered_train_X, y = filtered_train_E, weights = filtered_weights, 
                              degree = degree, nprune = nprune, glm = list(family = binomial),
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, ncross = ncross, minspan = minspan, endspan = endspan)
  }
  pred_onebin <- predict(fit.earth, newdata = newX, type = "response")
  pred <- 1 - (1 - pred_onebin)^num_intervals
  fit <- list(object = fit.earth, num_intervals = num_intervals)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("dIPCWSL.earth")
  return(out)
}


#' Prediction function for an dIPCWSL.earth object
#'
#' @param object Result object from dIPCWSL.earth
#' @param newdata Dataframe or matrix that will generate predictions
#' @param ... Any additional arguments (not used).
#' 
predict.dIPCWSL.earth <- function(object, newdata,...) {
  .SL.require('earth')
	pred_onebin <- predict(object$object, newdata = newdata, type = "response")
  num_intervals = object$num_intervals
  pred <- 1 - (1 - pred_onebin)^num_intervals
	return(pred)
}
