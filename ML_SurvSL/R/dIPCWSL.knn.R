#' Wrapper function for k-Nearest Neighbor models via yaImpute::ann().
#'
#' Note that the function only applies to binomial family, with outcomes bounded by [0, 1]
#'
#' @param Y Outcome variable
#' @param X Training dataframe
#' @param newX Test dataframe
#' @param time_point The cut-off time
#' @param family Binomial()
#' @param id Optional id to group observations from the same unit (not used currently).
#' @param num_intervals The number of equal time intervals into which the time_point should be divided.
#' @param k The number of nearest neighbors used to determine the prediction
#' @param ... Any remaining arguments, not used.
#' 
dIPCWSL.knn <- function(Y, X, newX, time_point, family, id, num_intervals = 1, k = 10, ...) {
  .SL.require('yaImpute')
  interval = seq(0, time_point, length.out = num_intervals + 1)
  expanded_train_data_with_weight = bin_combined_ipcw_separate(interval = interval, X = X, Y = Y,
                                                      time_point = time_point,return_type="dataset")
  varNames = colnames(X)
  train_weight_X = as.matrix(expanded_train_data_with_weight[varNames])
  newX = as.matrix(newX)
  fit.ann <- yaImpute::ann(ref = train_weight_X, target = newX, k = k, verbose = FALSE)
  neighbor_indices = as.data.frame(fit.ann$knnIndexDist)[, 1: k]
  neighbor_indices = as.matrix(neighbor_indices)
  pred_test_event_probability_single_bin = as.vector(apply(neighbor_indices, 1, function(neighbors) {
    E_values = expanded_train_data_with_weight$E[unlist(neighbors)]
    weights = expanded_train_data_with_weight$IPCW[unlist(neighbors)]
    above = sum(ifelse(weights == 0, 0, E_values * weights))
    below = sum(weights)
    c(est_test_event_single = ifelse(below == 0, 0, above/below))
  }))
  pred = 1 - (1 - pred_test_event_probability_single_bin)^num_intervals
  fit <- list(k = k, num_intervals = num_intervals)
  out <- list(pred = pred, fit=fit)
  class(out$fit) <- c("dIPCWSL.knn")
  return(out)
}


# will need original Y and X data for this
#' Prediction function for an dIPCWSL.knn object
#'
#' @param object Result object from dIPCWSL.knn
#' @param newdata Dataframe or matrix that will generate predictions.
#' @param X Dataframe of original X
#' @param Y Dataframe of original Y
#' @param ... Any additional arguments (not used).
#' 
predict.dIPCWSL.knn <- function(object, newdata, X, Y, ...){
  .SL.require('yaImpute')
  num_intervals = object$num_intervals
  interval = seq(0, time_point, length.out = num_intervals + 1)
  expanded_train_data_with_weight = bin_combined_ipcw_separate(interval = interval, X = X, Y = Y,
                                                               time_point = time_point,return_type="dataset")
  varNames = colnames(X)
  train_weight_X = as.matrix(expanded_train_data_with_weight[varNames])
  newdata = as.matrix(newdata)
  fit.ann <- yaImpute::ann(ref = train_weight_X, target = newdata, k = object$k, verbose = FALSE)
  neighbor_indices = as.data.frame(fit.ann$knnIndexDist)[, 1: object$k]
  neighbor_indices = as.matrix(neighbor_indices)
  pred_test_event_probability_single_bin = as.vector(apply(neighbor_indices, 1, function(neighbors) {
    E_values = expanded_train_data_with_weight$E[unlist(neighbors)]
    weights = expanded_train_data_with_weight$IPCW[unlist(neighbors)]
    above = sum(ifelse(weights == 0, 0, E_values * weights))
    below = sum(weights)
    c(est_test_event_single = ifelse(below == 0, 0, above/below))
  }))
  pred = 1 - (1 - pred_test_event_probability_single_bin)^num_intervals
  return(pred)
}
