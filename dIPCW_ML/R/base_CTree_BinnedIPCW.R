## ---------------------------------------------------------------------------------------------------------------
#' Build a weighted tree model and estimate event probabilities on X_validation using cross-validation results.
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
#' @param ... non-used parameters
#' 
base_CTree_BinnedIPCW <- function(X,
                                  Y,
                                  X_validation = NULL,
                                  ccp_alpha_values = NULL,
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
  if(is.null(ccp_alpha_values)){
    warning("No ccp values provided. Using grid: 2 * 10^seq(-1, -4, length = 40).")
    ccp_alpha_values = 2 * 10^seq(-1, -4, length = 40)
  }
  dt = data.frame(X, Y)
  var_name = colnames(X)
  returning = fit_classification_tree_binnedweights(X, 
                                                    Y, 
                                                    time_point = time_point,
                                                    measure = measure,
                                                    ccp_alpha_values = ccp_alpha_values,
                                                    range_intervals = range_intervals, 
                                                    k = k, 
                                                    foldid = foldid,
                                                    proxy_data = proxy_data,
                                                    ...)
  results = returning$results
  optimal_index_measure_row = get_optimal_row_indices_SingleMeasure(results, useMin = useMin, measure = measure)
  optimal_ind = optimal_index_measure_row[[1]]
  if(length(optimal_ind)!=1){
    print("CTreeError, optimal index should be length 1")
  }
  optimal_bins = results[optimal_ind, 1]
  optimal_ccp_alpha = results[optimal_ind, 2]
  interval = seq(0, time_point, length.out = optimal_bins + 1)
  expanded_train_data_with_weight = bin_combined_ipcw(interval = interval, data = dt, 
                                                      time_point = time_point, var_name = var_name, return_type="dataset")
  valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
  filtered_train_E = expanded_train_data_with_weight$E[valid_indices]
  filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
  formula = as.formula(paste("E ~", paste(var_name, collapse = " + ")))
  model_data = expanded_train_data_with_weight[valid_indices,]
  model_data$E = as.factor(model_data$E)
  tree_model = rpart(formula,
                     data = model_data,
                     weights = filtered_weights,
                     method = "class",
                     control = rpart.control(cp = optimal_ccp_alpha, minsplit = 5))
  
  if (!is.null(X_validation)) {
    pred_test_event_probability_single_bin = predict(tree_model, newdata = as.data.frame(X_validation), type = "prob")[,"1"]
    pred_test_surv_probability = (1 - pred_test_event_probability_single_bin)^optimal_bins
    pred_test_event_probability = 1 - pred_test_surv_probability
  } else {
    pred_test_event_probability = NULL
  }
  class_name = if (any(range_intervals != 1)) "Binned Tree" else "Single Tree"
  return(list(results = results, model = tree_model,
              pred = pred_test_event_probability, optbin = optimal_bins, class_name = class_name))
}
