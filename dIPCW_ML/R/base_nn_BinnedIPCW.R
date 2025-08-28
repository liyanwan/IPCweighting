#' Build a weighted Neural Network model and estimate event probabilities on X_validation using cross-validation results.
#'
#' @param X A dataframe of covariates.
#' @param Y A dataframe containing observed_time, sigma, and E.
#' @param X_validation A dataframe of covariates used for validation.
#' @param alpha A numeric vector specifying momentum rate for GD to evaluate.
#' @param lambda A numeric vector specifying L2 regularization parameter to evaluate.
#' @param lr_rate A numeric vector specifying learning rate parameter to evaluate.
#' @param time_point The specified cut-off time for evaluation.
#' @param measure The metric used for selecting the optimal tuning parameters.
#' @param range_intervals A numeric vector defining the number of intervals to evaluate during cross-validation.
#' @param k The number of folds used in cross-validation.
#' @param foldid A vector specifying the fold assignment for each observation.
#' @param useMin Logical; if TRUE, selects the tuning parameter based on the optimal measure value.
#' If FALSE, selects the parameter within one standard error of the optimal measure.
#' @param proxy_data A dataframe that may be used to compute IPCW in the test set.
#' @param sample_batch If TRUE, ensures that the same set of samples is selected in each batch across cross-validation folds.
#' @param batch_seed Used with set.seed(batch_seed) to guarantee reproducibility of batch sampling when sample_batch = TRUE.
#' 

base_deepnn_BinnedIPCW <- function(X, Y, X_validation = NULL,
                                   alpha, lambda, lr_rate,
                                   time_point,
                                   range_intervals,
                                   k = 5, 
                                   foldid = NULL,
                                   useMin = TRUE,
                                   proxy_data = NULL,
                                   model, epochs = 200, batch_size = 64, epsilon = 1e-3, verbose = FALSE,
                                   sample_batch = TRUE, batch_seed = 0,
                                   measure = c("C_index", "Log_Likelihood_Neg", "Brier_Score"),
                                   ...) {
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  dt = data.frame(X, Y)
  var_name = colnames(X)
  measure = match.arg(measure)
  cv = fit_deepnn_binnedweights(X = X, Y = Y, alpha = alpha, lambda = lambda, lr_rate = lr_rate, 
                                model = model, time_point = time_point, measure = measure,
                                k = k, foldid = foldid, proxy_data = proxy_data, epochs = epochs, 
                                batch_size = batch_size, epsilon = epsilon, verbose = verbose,
                                sample_batch = sample_batch, batch_seed = batch_seed,
                                range_intervals = range_intervals, ...)
  res = cv$results
  if(useMin){
    if(length(unique(cv$results[[measure]])) == 1){
      best = as.data.frame(t(vapply(cv$results,
                                     function(x) median(x, na.rm = TRUE),
                                     numeric(1))))
    } else {
      ord = if (identical(measure, "C_index")) order(-res[[measure]]) else order(res[[measure]])
      res = res[ord, , drop = FALSE]
      best = res[1, , drop = FALSE]
    }
  } else {
    warning("Invalid setting: useMin = FALSE. Only support useMin = TRUE.")
  }
  fit = NULL
  optimal_bins = best$num_bins
  interval = seq(0, time_point, length.out = optimal_bins + 1)
  expanded_train_data_with_weight = bin_combined_ipcw(interval = interval, data = dt, 
                                                      time_point = time_point, var_name = var_name, return_type="dataset")
  valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
  filtered_train_E = expanded_train_data_with_weight$E[valid_indices]
  filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
  formula = as.formula(paste("E ~", paste(var_name, collapse = " + ")))
  model_data = expanded_train_data_with_weight[valid_indices,]
  trainingX = as.matrix(model_data[, var_name, drop = FALSE])
  if(sample_batch) set.seed(batch_seed)
  fit = dnn::deepGlm(filtered_train_E ~ trainingX,
                      model = model, family = "binomial", weights = filtered_weights,
                      epochs = epochs, lr_rate = best$lr_rate, alpha = best$alpha, lambda = best$lambda,
                      batch_size = batch_size, verbose = verbose)
  if (!is.null(X_validation)) {
    est_test_event_single = as.numeric(predict(fit$model, x = as.matrix(X_validation)))
    est_test_surv_probability = (1 - est_test_event_single)^num_intervals
    est_test_event_probability = 1 - est_test_surv_probability
  } else {
    est_test_event_probability = NULL
  }
  class_name = paste0(if (any(range_intervals != 1)) "Binned " else "Single ", "NNet")
  return(list(results = res, model = fit, 
              pred = est_test_event_probability, optbin = optimal_bins, class_name = class_name))
}
