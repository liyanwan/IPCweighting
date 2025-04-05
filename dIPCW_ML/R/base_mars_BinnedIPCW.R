#' Build a weighted MARS model and estimate event probabilities on X_validation using cross-validation results.
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
#' @param var_threshold The number of covariates selected in earth modelling. if NULL, select all covariates
#' @param ... Any additional arguments are passed through to earth.
#' 
base_mars_BinnedIPCW <- function(X,
                                 Y,
                                 X_validation = NULL,
                                 nprune_values = NULL,
                                 degree_values = NULL,
                                 time_point,
                                 range_intervals,
                                 k = 5, 
                                 foldid = NULL,
                                 useMin = FALSE,
                                 proxy_data = NULL,
                                 measure = "C_index",
                                 # var_threshold = NULL,
                                 ...){
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  if(is.null(nprune_values)){
    warning("No nprune values provided. Using the default value NULL.")
  }
  if(is.null(degree_values)){
    warning("No degree values provided. Using value 2")
    degree_values = 2
  }
  if(is.null(X_validation)){
    X_validation = X
  }
  dt = data.frame(X, Y)
  var_name = colnames(X)
  returning = fit_mars_binnedweights(X,
                                   Y,
                                   time_point,
                                   nprune_values = nprune_values,
                                   degree_values = degree_values,
                                   range_intervals = range_intervals,
                                   k = k, foldid = foldid,
                                   proxy_data = proxy_data,
                                   measure = measure,
                                   # var_threshold = var_threshold,
                                   ...)
  results = returning$results
    optimal_index_measure_row = get_optimal_row_indices_SingleMeasure(results, useMin = useMin, measure = measure)
    optimal_ind = optimal_index_measure_row[[1]]
    optimal_bins = results[optimal_ind, 1]
    optimal_nprune = results[optimal_ind, 2]
    optimal_degree = results[optimal_ind, 3]
    interval = seq(0, time_point, length.out = optimal_bins + 1)
    expanded_train_data_with_weight = bin_combined_ipcw(interval,dt, time_point, var_name, return_type="dataset")
    valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
    filtered_train_E = expanded_train_data_with_weight$E[valid_indices]
    filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
    rss_vals = numeric(length(var_name))
    model_data = expanded_train_data_with_weight[valid_indices,]
    # for (vi in seq_along(var_name)) {
    #   variable_name = var_name[vi]
    #   form = as.formula(paste("E ~", variable_name))
    #   fit_uni = earth(form,
    #                   data = model_data,
    #                   glm = list(family = binomial),
    #                   weights = filtered_weights)
    #   rss_vals[vi] = fit_uni$rss
    # }
    # sig_var_name = var_name[rank(rss_vals) <= var_threshold]
    sig_var_name = var_name
    formula = as.formula(paste("E ~", paste(sig_var_name, collapse = " + ")))
    earth_model = earth(formula,
                        data = model_data,
                        weights = filtered_weights,
                        glm = list(family = binomial),
                        nprune = optimal_nprune,
                        degree = optimal_degree)
    pred_test_event_probability_single_bin = predict(earth_model, newdata = X_validation, type = "response")
    pred_test_surv_probability = (1 - pred_test_event_probability_single_bin)^optimal_bins
    pred_test_event_probability = 1 - pred_test_surv_probability
    
  class_name = if (any(range_intervals != 1)) "Binned MARS" else "Single MARS"
  return(list(results = results, model = earth_model, 
              pred = pred_test_event_probability, optbin = optimal_bins, class_name = class_name))
}
