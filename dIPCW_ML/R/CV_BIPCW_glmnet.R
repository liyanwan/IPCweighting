#' Compute Cross-Validation Results Using an IPC weighted elastic-net Model with discrete-time modelling
#'
#' @param X A dataframe of covariates.
#' @param Y A dataframe containing observed_time, sigma, and E.
#' @param time_point The specified cut-off time for evaluation.
#' @param measure The metric used for selecting the optimal tuning parameters.
#' @param alpha Elastic net mixing parameter, range [0, 1]. 0 = ridge regression
#'   and 1 = lasso.
#' @param lambda_grid A vector of lambda to check.
#' @param range_intervals A numeric vector specifying the number of intervals to evaluate during cross-validation.
#' @param k The number of cross-validation folds.
#' @param foldid A vector indicating the fold assignment for each observation.
#' @param proxy_data A dataframe that may be used to compute IPCW in the test set.
#' @param ... Any additional arguments are passed through to glmnet.

fit_glmnet_binnedweights <- function(X, 
                                    Y, 
                                    time_point, 
                                    measure,
                                    alpha_num, 
                                    lambda_values, 
                                    range_intervals, 
                                    k = 5, 
                                    foldid = NULL,
                                    proxy_data = NULL, 
                                    ...) {
  dt = cbind(X, Y)
  var_name = colnames(X)
  n = nrow(dt)
  nlambda = length(lambda_values)
  cv_results = expand.grid(num_bins = range_intervals,
                          lambda = lambda_values)
  cv_results[[measure]] = NA
  if(is.null(foldid)){
    foldid = sample(rep(1:k, length.out = n))
  }
  if(is.list(foldid)){
    stop("foldid should be a vector of integers")
  }
  k = max(foldid)
  for (i in range_intervals) {
    num_intervals <- i
    fold_performance_measure = array(NA, dim = c(k, length(lambda_values)))
    interval = seq(0, time_point, length.out = num_intervals + 1)
    for (f in 1:k) {
      train_data = dt[foldid != f, ]
      test_data = dt[foldid == f, ]
      test_X = test_data[var_name]
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
      expanded_train_data_with_weight = bin_combined_ipcw(interval = interval, data = train_data,
                                         var_name = var_name, time_point = time_point, return_type = "dataset")
      train_weight_X = expanded_train_data_with_weight[var_name]
      valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
      filtered_train_X = train_weight_X[valid_indices, ]
      filtered_train_E = expanded_train_data_with_weight$E[valid_indices]
      filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
      penalty_model <- tryCatch({
        glmnet(
          as.matrix(filtered_train_X),
          filtered_train_E,
          family = "binomial",
          weights = filtered_weights,
          alpha = alpha_num,
          lambda = lambda_values,
          standardize = TRUE,
          ...
        )
      }, error = function(e) {
        message("Error in fitting glmnet for lambda = ", lambda_values, ": ", e$message)
        return(NULL)
      })
      
      if (!is.null(penalty_model)) {
        est_test_surv_probability = (1 - predict(penalty_model,
                                                 newx = as.matrix(test_X),
                                                 type = "response"))^num_intervals
        est_test_event_probability = 1 - est_test_surv_probability
        
        for (j in 1:ncol(est_test_event_probability)) {
          fold_performance_measure[f, j] <- switch(measure,
                                                   "C_index" = concordance.index(x = est_test_event_probability[, j],
                                                                                 surv.time = test_data$observed_time,
                                                                                 surv.event = test_data$sigma)$c.index,
                                                   "Log_Likelihood_Neg" = weighted_loglikelihood(IPCW = test_IPCW, 
                                                                                                 status = test_data$E, 
                                                                                                 event_prob = est_test_event_probability[, j]),
                                                   "Brier_Score" = Weighted_Brier_Score(event_prob = est_test_event_probability[, j], 
                                                                                        status = test_data$E, 
                                                                                        IPCW = test_IPCW),
                                                   stop("Invalid measure. Please use 'C_index', 'Log_Likelihood_Neg', or 'Brier_Score'.")
                                                   )
          if (measure == "C_index" && is.na(fold_performance_measure[f, j])) 
            fold_performance_measure[f, j] = 0.5
        }
        npreds = ncol(est_test_event_probability)
        if(npreds < nlambda){
          fold_performance_measure[f,seq(from=npreds,to=nlambda)]=fold_performance_measure[f,npreds]
        }
      }
    }
    fold_sizes = table(foldid)
    for (j in 1:length(lambda_values)) {
      valid_indices <- !is.na(fold_performance_measure[, j])
      fold_weights = fold_sizes[valid_indices]/sum(fold_sizes[valid_indices])
      cv_results[[measure]][cv_results$num_bins == num_intervals & cv_results$lambda == lambda_values[j]] = sum(fold_performance_measure[valid_indices, j] * fold_weights)
      }
    }
  return(list(optimal_results = cv_results))
}
