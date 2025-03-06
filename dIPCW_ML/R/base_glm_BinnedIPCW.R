## ---------------------------------------------------------------------------------------------------------------
#' Build a weighted GLM model and estimate event probabilities on X_validation using cross-validation results.
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

base_glm_BinnedIPCW <- function(X,
                                Y,
                                X_validation = NULL,
                                time_point,
                                range_intervals,
                                k = 5,
                                foldid = NULL,
                                useMin = FALSE,
                                proxy_data = NULL,
                                measure = "C_index"){
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  if(is.null(X_validation)){
    X_validation = X
  }
  dt = data.frame(X, Y)
  var_name = colnames(X)
  if(length(range_intervals) != 1){
    returning = fit_glm_binnedweights(X = X,
                                       Y = Y,
                                       k = k,
                                       foldid = foldid,
                                       measure = measure,
                                       time_point = time_point,
                                       range_intervals = range_intervals,
                                       proxy_data = proxy_data)
    results = returning$results
    optimal_index_all_row = get_optimal_row_indices_SingleMeasure(results, useMin = useMin, measure = measure)
  } else {
    optimal_index_all_row = list(1)
    results = data.frame(num_bins = range_intervals)
  }
  ind = optimal_index_all_row[[1]]
  if(length(ind) != 1){
    stop("GLM error: optimal index should be length 1")
  }
  optimal_bins = results[ind, 1]
  interval_ind = seq(0, time_point, length.out = optimal_bins + 1)
  easy_model_ind = bin_combined_ipcw(interval = interval_ind, data = dt, time_point = time_point, var_name = var_name)
  if (!is.null(X_validation)) {
    X_validation = as.data.frame(X_validation)
    pred_test_event_probability_single_bin = predict(easy_model_ind,
                                                     newdata = X_validation,
                                                     type = "response")
    pred_test_surv_probability = (1 - pred_test_event_probability_single_bin)^optimal_bins
    pred_test_event_probability = 1 - pred_test_surv_probability
  } else {
    pred_test_event_probability = NULL
  }
  class_name = if (any(range_intervals != 1)) "Binned GLM" else "Single GLM"
  return(list(results = results, pred = pred_test_event_probability,
              model = easy_model_ind, optbin = optimal_bins, class_name = class_name))
}
