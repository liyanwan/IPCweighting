## ---------------------------------------------------------------------------------------------------------------
#' Compute Cross-Validation Results Using an IPC weighted GLM Model with discrete-time modelling
#'
#' @param X A dataframe of covariates.
#' @param Y A dataframe containing observed_time, sigma, and E.
#' @param time_point The specified cut-off time for evaluation.
#' @param measure The metric used for selecting the optimal tuning parameters.
#' @param range_intervals A numeric vector specifying the number of intervals to evaluate during cross-validation.
#' @param k The number of cross-validation folds.
#' @param foldid A vector indicating the fold assignment for each observation.
#' @param proxy_data A dataframe that may be used to compute IPCW in the test set.

fit_glm_binnedweights <- function(X, 
                          Y, 
                          time_point, 
                          measure,
                          range_intervals, 
                          k = 5, 
                          foldid = NULL,
                          proxy_data = NULL){
  dt = cbind(X, Y)
  var_name = colnames(X)
  results = data.frame(intervals = sort(range_intervals))
  results[[measure]] = NA
  if(is.null(foldid)){
    foldid = sample(rep(1:k, length.out = nrow(dt)))
  }
  k = max(foldid)
  fold_sizes = table(foldid)
  fold_weights = fold_sizes/sum(fold_sizes)
  
  for (num_intervals in sort(range_intervals)) {
    fold_performance_measure = numeric(k)
    for (f in 1:k) {
      train_data = dt[foldid != f, ]
      test_data = dt[foldid == f, ]
      train_Ghat = Get_Ghat(observed_time = train_data$observed_time, sigma = train_data$sigma)
      # Compute IPC weights of test data using G from train data
      if (!is.null(proxy_data)) {
        test_proxy = rbind(test_data, proxy_data)
        totaltest_G = Ghat_newtime(step_function = train_Ghat, 
                                   new_observed_time = test_proxy$observed_time, 
                                   time_point = time_point)
        total_IPCW = ifelse(is.na(test_data$E), 0, 1/totaltest_G)
        test_IPCW = total_IPCW[1:nrow(test_data)]
      } else {
        test_G = Ghat_newtime(step_function = train_Ghat,
                              new_observed_time = test_data$observed_time,
                              time_point = time_point)
        test_IPCW = ifelse(is.na(test_data$E), 0, 1/test_G)
      }
      interval = seq(0, time_point, length.out = num_intervals + 1)
      model = bin_combined_ipcw(interval = interval, data = train_data,
                                var_name = var_name, time_point = time_point, return_type = "model")
      pred_combine_event_prob = predict(model, newdata = test_data[var_name], type = "response")
      pred_combine_surv_prob = 1 - pred_combine_event_prob
      pred_risk_prob = 1-(pred_combine_surv_prob)^num_intervals
      fold_performance_measure[f] <- switch(measure,
                                            "C_index" = concordance.index(x = pred_risk_prob,
                                                                          surv.time = test_data$observed_time,
                                                                          surv.event = test_data$sigma)$c.index,
                                            "Log_Likelihood_Neg" = weighted_loglikelihood(IPCW = test_IPCW, 
                                                                                          status = test_data$E, 
                                                                                          event_prob = pred_risk_prob),
                                            "Brier_Score" = Weighted_Brier_Score(event_prob = pred_risk_prob, 
                                                                                 status = test_data$E, 
                                                                                 IPCW = test_IPCW),
                                            stop("Invalid measurement type. Please try another one.")
                                            )
    }
    results[[measure]][num_intervals] = sum(fold_performance_measure * fold_weights)
  }
  return(list(results = results, measure = measure))
}
