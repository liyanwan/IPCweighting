## ---------------------------------------------------------------------------------------------------------------
#' Build a weighted elastic net model and estimate event probabilities on X_validation using cross-validation results.
#'
#' @param X A dataframe of covariates.
#' @param Y A dataframe containing observed_time, sigma, and E.
#' @param X_validation A dataframe of covariates used for validation.
#' @param time_point The specified cut-off time for evaluation.
#' @param measure The metric used for selecting the optimal tuning parameters.
#' @param range_intervals A numeric vector defining the number of intervals to evaluate during cross-validation.
#' @param k The number of folds used in cross-validation.
#' @param foldid A vector specifying the fold assignment for each observation.
#' @param useMin Logical; if TRUE, selects the tuning parameter based on the optimal measure value.
#' If FALSE, selects the parameter within one standard error of the optimal measure.
#' @param proxy_data A dataframe that may be used to compute IPCW in the test set.
#' @param ... Any additional arguments are passed through to glmnet.

base_glmnet_BinnedIPCW <- function(X, 
                                   Y, 
                                   X_validation = NULL,
                                   lambda_values = NULL,
                                   alpha_num = 1,
                                   time_point,
                                   range_intervals,
                                   k = 5, 
                                   foldid = NULL,
                                   useMin = FALSE,
                                   proxy_data = NULL,
                                   measure = "C_index",
                                   ...){
  
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  if(is.null(lambda_values)){
    warning("No lambda values provided. Using grid: 10^seq(-1, -5, length = 50).")
    lambda_values = 10^seq(-1, -5, length = 50)
  }
  dt = data.frame(X, Y)
  var_name = colnames(X)
  results = fit_glmnet_binnedweights(X,
                                    Y,
                                    time_point,
                                    alpha_num = alpha_num,
                                    lambda_values = lambda_values,
                                    range_intervals = range_intervals,
                                    k = k, 
                                    foldid = foldid,
                                    proxy_data = proxy_data,
                                    measure = measure,
                                    ...) # glmnet parameters
  optimal_results = results$optimal_results
  optimal_index_measure_rows = get_optimal_row_indices_SingleMeasure(optimal_results, useMin = useMin, measure = measure)
  optimal_row = optimal_index_measure_rows[[1]]
  if(length(optimal_row) != 1){
    stop("GLMNET error: optimal index should be length 1")
  }
  optimal_bins = optimal_results[optimal_row, 1]
  optimal_lambda = optimal_results[optimal_row, 2]
  interval = seq(0, time_point, length.out = optimal_bins + 1)
  expanded_train_data_with_weight = bin_combined_ipcw(interval = interval, data = dt, time_point = time_point, 
                                                      var_name = var_name, return_type="dataset")
  weight_X = expanded_train_data_with_weight[var_name]
  valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
  filtered_X = weight_X[valid_indices, ]
  filtered_E = expanded_train_data_with_weight$E[valid_indices]
  filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
  lambda_penalty_model <- glmnet(filtered_X, 
                                 filtered_E, 
                                 family = "binomial", 
                                 weights = filtered_weights, 
                                 alpha = alpha_num, 
                                 lambda = lambda_values,...)
    if (!is.null(X_validation)) {
      pred_test_event_probability_single_bin = predict(lambda_penalty_model, 
                                                       newx = as.matrix(X_validation), 
                                                       s = optimal_lambda,
                                                       type = "response")
      pred_test_surv_probability = (1 - pred_test_event_probability_single_bin)^optimal_bins
      pred_test_event_probability = 1 - pred_test_surv_probability
    } else {
      pred_test_event_probability = NULL
    }
  penalized_type = switch(
    as.character(alpha_num),
    "0" = "Ridge",
    "1" = "Lasso",
    "Elastic"
  )
  class_name = paste0(if (any(range_intervals != 1)) "Binned " else "Single ", penalized_type)
  return(list(results = results, model = lambda_penalty_model, 
              pred = pred_test_event_probability, optbin = optimal_bins, class_name = class_name))
}
