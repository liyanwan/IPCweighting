#' Wrapper function for generalized linear models via glm().
#'
#' Note that the function only applies to binomial family, with outcomes bounded by [0, 1]
#'
#' @param Y Outcome variable
#' @param X Training dataframe
#' @param newX Test dataframe
#' @param time_point The cut-off time
#' @param family Binomial()
#' @param num_intervals The number of equal time intervals into which the time_point should be divided.
#' @param ... Any remaining arguments, not used.
#'
#' @references
#' Polley et al., 2011.

dIPCWSL.glm <- function(Y, X, newX, time_point, family, num_intervals = 1, ...) {

  # X must be a dataframe, not a matrix.
  if (is.matrix(X)) {
    X <- as.data.frame(X)
  }
  interval <- seq(0, time_point, length.out = num_intervals + 1)
  fit.glm <- bin_combined_ipcw_separate(interval = interval, X = X, Y = Y,
                               time_point = time_point, return_type="model")
  # newX must be a dataframe, not a matrix.
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred_onebin <- predict(fit.glm, newdata = newX, type = "response")
  pred <- 1 - (1 - pred_onebin)^num_intervals
  fit <- list(object = fit.glm, num_intervals = num_intervals)
  class(fit) <- "dIPCWSL.glm"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction function for dIPCWSL.glm
#'
#' @param object dIPCWSL.glm object
#' @param newdata Dataframe to generate predictions
#' @param ... Unused additional arguments

predict.dIPCWSL.glm <- function(object, newdata, ...) {
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  pred_onebin <- predict(object = object$object, newdata = newdata, type = "response")
  num_intervals = object$num_intervals
  pred <- 1 - (1 - pred_onebin)^num_intervals
  return(pred)
}
