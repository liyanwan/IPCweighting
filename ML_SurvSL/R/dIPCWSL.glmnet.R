#' Wrapper function using elastic net. Alpha = 0 corresponds to ridge, Alpha = 1 corresponds to lasso
#'
#' @param Y Outcome variable
#' @param X Covariate dataframe
#' @param newX Dataframe to predict the outcome
#' @param time_point cut-off time
#' @param family "binomial" for binary classification.
#' @param id Optional id to group observations from the same unit (not used
#'   currently).
#' @param num_intervals The number of equal time intervals into which the time_point should be divided.
#' @param alpha Elastic net mixing parameter, range [0, 1]. 0 = ridge regression
#'   and 1 = lasso.
#' @param lambda_grid A vector of lambda to check.
#' @param ... Any additional arguments are passed through to glmnet.
#' @seealso \code{\link{SL.glmnet}} \code{\link[glmnet]{predict.dIPCWSL.glmnet}}
#'   \code{\link[glmnet]{glmnet}}

dIPCWSL.glmnet <- function(Y, X, newX, time_point, family, id,
                           num_intervals = 1, alpha = 1, lambda_grid = NULL, ...) {
  .SL.require('glmnet')
  # X must be a matrix, should we use model.matrix or as.matrix
  if (!is.matrix(X)) {
    X <- model.matrix(~ -1 + ., X)
    newX <- model.matrix(~ -1 + ., newX)
  }
  interval <- seq(0, time_point, length.out = num_intervals + 1)
  expanded_train_data_with_weight <- bin_combined_ipcw_separate(interval = interval, X = X, Y = Y,
                                                       time_point = time_point,  return_type = "dataset")
  varNames = colnames(X)
  train_weight_X <- expanded_train_data_with_weight[varNames]
  valid_indices <- which(expanded_train_data_with_weight$IPCW > 0)
  filtered_train_X <- train_weight_X[valid_indices, ]
  filtered_train_E <- expanded_train_data_with_weight$E[valid_indices]
  filtered_weights <- expanded_train_data_with_weight$IPCW[valid_indices]

  fit.glmnet <- tryCatch({
    glmnet::glmnet(
      x = as.matrix(filtered_train_X),
      y = filtered_train_E,
      weights = filtered_weights,
      lambda = lambda_grid,
      family = family$family,
      alpha = alpha,
      ...
    )
  }, error = function(e) {
    message("Error in glmnet: ", e$message)
    return(NULL)
  })

  # If we predict with the cv.glmnet object we can specify lambda using a string.
  if(is.null(fit.glmnet)){
    pred_onebin = rep(NA, nrow(newX))
  } else{
    pred_onebin = predict(fit.glmnet, newx = newX, type = "response")
  }
  pred <- 1 - (1 - pred_onebin)^num_intervals
  fit <- list(object = fit.glmnet, num_intervals = num_intervals)
  class(fit) <- "dIPCWSL.glmnet"

  out <- list(pred = pred, fit = fit)
  return(out)
}

#' Prediction function for an dIPCWSL.glmnet object
#'
#' @param object Result object from dIPCWSL.glmnet
#' @param newdata Dataframe or matrix that will generate predictions.
#' @param remove_extra_cols Remove any extra columns in the new data that were
#'   not part of the original model.
#' @param add_missing_cols Add any columns from original data that do not exist
#'   in the new data, and set values to 0.
#' @param ... Any additional arguments (not used).

predict.dIPCWSL.glmnet <- function(object, newdata,
                              remove_extra_cols = T,
                              add_missing_cols = T,
                              ...) {
  .SL.require('glmnet')

  if (!is.matrix(newdata)) {
    newdata <- model.matrix(~ -1 + ., newdata)
  }
  original_cols = rownames(object$object$glmnet.fit$beta)
  num_intervals = object$num_intervals

  # Remove any columns in newdata that were not present in original data.
  if (remove_extra_cols) {
    extra_cols = setdiff(colnames(newdata), original_cols)
    if (length(extra_cols) > 0) {
      warning(paste("Removing extra columns in prediction data:",
                     paste(extra_cols, collapse = ", ")))

      newdata = newdata[, !colnames(newdata) %in% extra_cols, drop = FALSE]
    }
  }

  # Add any columns in original data that are not present in new data.
  if (add_missing_cols) {
    missing_cols = setdiff(original_cols, colnames(newdata))
    if (length(missing_cols) > 0) {
      warning(paste("Adding missing columns in prediction data:",
                     paste(missing_cols, collapse = ", ")))

      new_cols = matrix(0, nrow = nrow(newdata), ncol = length(missing_cols))
      colnames(new_cols) = missing_cols
      newdata = cbind(newdata, new_cols)

      # Sort columns in the correct order so that matrix multiplication is correct.
      newdata = newdata[, original_cols, drop = FALSE]
    }
  }
  pred_onebin <- predict(object$object, newx = newdata, type = "response")
  pred <- 1 - (1 - pred_onebin)^num_intervals

  return(pred)
}
