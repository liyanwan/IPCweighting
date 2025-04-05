## ---------------------------------------------------------------------------------------------------------------
generate_Y_continuous <- function(X, ind, beta) {
  if (length(ind) != length(beta)) {
    stop("Length of 'ind' and 'beta' must be the same.")
  }
  Y_continuous = X[, ind] %*% beta
  return(Y_continuous)
}


## ---------------------------------------------------------------------------------------------------------------
predict_nodes <- function (object, newdata, na.action = na.pass) {
  where <-
    if (missing(newdata)) 
      object$where
  else {
    if (is.null(attr(newdata, "terms"))) {
      Terms <- delete.response(object$terms)
      newdata <- model.frame(Terms, newdata, na.action = na.action, 
                             xlev = attr(object, "xlevels"))
      if (!is.null(cl <- attr(Terms, "dataClasses"))) 
        .checkMFClasses(cl, newdata, TRUE)
    }
    rpart:::pred.rpart(object, rpart:::rpart.matrix(newdata))
  }
  as.integer(row.names(object$frame))[where]
}


## ---------------------------------------------------------------------------------------------------------------
#' Compute Cross-Validation Results Using an IPC weighted Tree Model with discrete-time modelling
#'
#' @param X A dataframe of covariates.
#' @param Y A dataframe containing observed_time, sigma, and E.
#' @param time_point The specified cut-off time for evaluation.
#' @param measure The metric used for selecting the optimal tuning parameters.
#' @param ccp_alpha_values A numeric vector specifying cp parameter to evaluate.
#' @param range_intervals A numeric vector specifying the number of intervals to evaluate during cross-validation.
#' @param k The number of cross-validation folds.
#' @param foldid A vector indicating the fold assignment for each observation.
#' @param proxy_data A dataframe that may be used to compute IPCW in the test set.
#' 
fit_classification_tree_binnedweights <- function(X, 
                                                  Y, 
                                                  time_point,
                                                  measure,
                                                  ccp_alpha_values,
                                                  range_intervals, 
                                                  k = 5, 
                                                  foldid = NULL,
                                                  proxy_data = NULL,
                                                  ...) { 
  data = cbind(X, Y)
  var_name = colnames(X)
  cv_results <- expand.grid(
    num_bins = range_intervals,
    ccp_alpha_values = ccp_alpha_values
  )
  cv_results[[measure]] = NA
  if(is.null(foldid)){
    foldid = sample(rep(1:k, length.out = nrow(X)))
  }
  k = max(foldid)
  fold_sizes = table(foldid)
  fold_weights = fold_sizes/sum(fold_sizes)
  
  for (i in 1:nrow(cv_results)) {
    num_intervals = cv_results$num_bins[i]
    interval = seq(0, time_point, length.out = num_intervals + 1)
    ccp_alpha = cv_results$ccp_alpha_values[i]
    fold_performance_measure = numeric(k)
    
    for (f in 1:k) {
      train_data = data[foldid != f, ]
      test_data = data[foldid == f, ]
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
      expanded_train_data_with_weight = bin_combined_ipcw(interval = interval,data = train_data,
                                                          var_name = var_name, time_point = time_point, 
                                                          return_type="dataset")
      valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
      filtered_train_E = expanded_train_data_with_weight$E[valid_indices]
      filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
      formula = as.formula(paste("E ~", paste(var_name, collapse = " + ")))
      model_data = expanded_train_data_with_weight[valid_indices,]
      tree_model = rpart(formula,
                         data = model_data,
                         weights = filtered_weights,
                         method = "class",
                         control = rpart.control(cp = ccp_alpha))
      train_node_indices = predict_nodes(tree_model, model_data)
      test_node_indices = predict_nodes(tree_model, test_data)
      # est_test_label = predict(tree_model, newdata = test_data, type = "class")
      node_probabilities = tapply(X = filtered_weights * filtered_train_E,
                                  INDEX = train_node_indices,
                                  FUN = sum) / tapply(filtered_weights, INDEX=train_node_indices, FUN = sum)
      node_probabilities = setNames(node_probabilities, names(node_probabilities))
      est_test_event_single = node_probabilities[as.character(test_node_indices)]
      est_test_surv_probability = (1-est_test_event_single)^num_intervals
      est_test_event_probability = 1-est_test_surv_probability
      fold_performance_measure[f] <- switch(measure,
                                            "C_index" = concordance.index(x = est_test_event_probability,
                                                                          surv.time = test_data$observed_time,
                                                                          surv.event = test_data$sigma)$c.index,
                                            "Log_Likelihood_Neg" = weighted_loglikelihood(IPCW = test_IPCW, 
                                                                                          status = test_data$E, 
                                                                                          event_prob = est_test_event_probability),
                                            "Brier_Score" = Weighted_Brier_Score(event_prob = est_test_event_probability, 
                                                                                 status = test_data$E, 
                                                                                 IPCW = test_IPCW),
                                            stop("Invalid measure. Please use 'C_index', 'Log_Likelihood_Neg', or 'Brier_Score'.")
                                            )
    }
    if(measure == "C_index"){
      valid_indices_c <- !is.na(fold_performance_measure)
      fold_weights_c = fold_sizes[valid_indices_c]/sum(fold_sizes[valid_indices_c])
      cv_results[[measure]][i] = sum(fold_performance_measure[valid_indices_c]*fold_weights_c)
    }
    cv_results[[measure]][i] = sum(fold_performance_measure * fold_weights)
  }
  return(list(results = cv_results, measure = measure))
}
