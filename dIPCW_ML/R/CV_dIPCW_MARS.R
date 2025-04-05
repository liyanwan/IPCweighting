#' Compute Cross-Validation Results Using a IPC weighted mars Model with discrete-time modelling
#'
#' @param X A dataframe of covariates.
#' @param Y A dataframe containing observed_time, sigma, and E.
#' @param time_point The specified cut-off time for evaluation.
#' @param measure The metric used for selecting the optimal tuning parameters.
#' @param nprune_values A numeric vector specifying the nprune parameter values to evaluate in cross-validation
#' @param degree_values A numeric vector specifying the degree parameter values to evaluate in cross-validation
#' @param range_intervals A numeric vector specifying the number of intervals to evaluate in cross-validation.
#' @param k The number of cross-validation folds.
#' @param foldid A vector indicating the fold assignment for each observation.
#' @param proxy_data A dataframe that may be used to compute IPCW in the test set.
#' @param var_threshold The number of covariates selected in earth modelling. if NULL, select all covariates
#' @param ... Any additional arguments are passed through to earth.
#' 
# fit_mars_binnedweights <- function(X,
#                              Y,
#                              time_point,
#                              measure,
#                              nprune_values,
#                              degree_values,
#                              range_intervals,
#                              k = 5,
#                              foldid = NULL,
#                              proxy_data = NULL,
#                              var_threshold = NULL,
#                              ...) {
#   data = cbind(X, Y)
#   var_name = colnames(X)
#   if(is.null(nprune_values)){
#     nprune_values = NA
#   }
#   cv_results = expand.grid(
#     num_bins = range_intervals,
#     nprune_values = nprune_values,
#     degree_values = degree_values
#   )
#   cv_results[[measure]] = NA
#   if(is.null(foldid)){
#     foldid = sample(rep(1:k, length.out = nrow(data)))
#   }
#   if(is.null(var_threshold)){
#     var_threshold = ncol(X)
#   }
#   k = max(foldid)
#   fold_sizes = table(foldid)
#   fold_weights = fold_sizes/sum(fold_sizes)
# 
#   for (i in 1:(nrow(cv_results))) {
#     fold_performance_measure = numeric(k)
#     num_intervals = cv_results$num_bins[i]
#     interval = seq(0, time_point, length.out = num_intervals + 1)
#     degree = cv_results$degree_values[i]
#     nprune = cv_results$nprune_values[i]
#     if(is.na(nprune)){
#       nprune = NULL
#     }
#     for (j in 1:k) {
#       train_data = data[foldid != j, ]
#       test_data = data[foldid == j, ]
#       train_Ghat = Get_Ghat(observed_time = train_data$observed_time, sigma = train_data$sigma)
#       # Compute IPC weights of test data using G from train data
#       if (!is.null(proxy_data)) {
#         test_proxy = rbind(test_data, proxy_data)
#         totaltest_G = Ghat_newtime(step_function = train_Ghat,
#                                    new_observed_time = test_proxy$observed_time,
#                                    time_point = time_point)
#         total_IPCW = ifelse(test_data$sigma == 0 & test_data$observed_time<=time_point, 0, 1/totaltest_G)
#         test_IPCW = total_IPCW[1:nrow(test_data)]
#       } else {
#         test_G = Ghat_newtime(step_function = train_Ghat,
#                               new_observed_time = test_data$observed_time,
#                               time_point = time_point)
#         test_IPCW = ifelse(test_data$sigma == 0 & test_data$observed_time<=time_point, 0, 1/test_G)
#       }
#       # Get the expanded train dataset using combined IPCW
#       expanded_train_data_with_weight = bin_combined_ipcw(interval = interval,data = train_data,
#                                                           var_name = var_name, time_point = time_point,
#                                                           return_type="dataset")
#       valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
#       filtered_train_E = expanded_train_data_with_weight$E[valid_indices]
#       filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
#       # Fit MARS with BINNED IPCW weights
#       rss_vals = numeric(length(var_name))
#       model_data = expanded_train_data_with_weight[valid_indices,]
#       for (vi in seq_along(var_name)) {
#         variable_name = var_name[vi]
#         form = as.formula(paste("E ~", variable_name))
#         fit_uni = earth(form,
#                         data = model_data,
#                         glm = list(family = binomial),
#                         weights = filtered_weights)
#         rss_vals[vi] = fit_uni$rss
#       }
#       sig_var_name = var_name[rank(rss_vals) <= var_threshold]
#       formula = as.formula(paste("E ~", paste(sig_var_name, collapse = " + ")))
#       earth_model = earth(formula,
#                           data = model_data,
#                           weights = filtered_weights,
#                           glm = list(family = binomial),
#                           nprune = nprune,
#                           degree = degree,
#                           ...)
#       est_test_event_single = predict(earth_model, newdata = test_data, type = "response")
#       est_test_surv_probability = (1-est_test_event_single)^num_intervals
#       est_test_event_probability = 1-est_test_surv_probability
#       fold_performance_measure[j] <- switch(measure,
#                                             "C_index" = concordance.index(x = est_test_event_probability,
#                                                                           surv.time = test_data$observed_time,
#                                                                           surv.event = test_data$sigma)$c.index,
#                                             "Log_Likelihood_Neg" = weighted_loglikelihood(IPCW = test_IPCW,
#                                                                                           status = test_data$E,
#                                                                                           est_test_event_probability),
#                                             "Brier_Score" = Weighted_Brier_Score(est_test_event_probability,
#                                                                                  status = test_data$E,
#                                                                                  IPCW = test_IPCW),
#                                             stop("Invalid measure. Please use 'C_index', 'Log_Likelihood_Neg', or 'Brier_Score'.")
#                                             )
#     }
#     cv_results[[measure]][i] = sum(fold_performance_measure * fold_weights)
#   }
#   return(list(results = cv_results))
# }

fit_mars_binnedweights <- function(X,
                             Y,
                             time_point,
                             measure,
                             nprune_values,
                             degree_values,
                             range_intervals,
                             k = 5,
                             foldid = NULL,
                             proxy_data = NULL,
                             # var_threshold = NULL,
                             ...) {
  data = cbind(X, Y)
  var_name = colnames(X)
  if(is.null(nprune_values)){
    nprune_values = NA
  }
  cv_results = expand.grid(
    num_bins = range_intervals,
    nprune_values = nprune_values,
    degree_values = degree_values
  )
  cv_results[[measure]] = NA
  cross_validation_grid = expand.grid(
    num_bins = range_intervals,
    degree_values = degree_values
  )
  if(is.null(foldid)){
    foldid = sample(rep(1:k, length.out = nrow(data)))
  }
  if(is.null(var_threshold)){
    var_threshold = ncol(X)
  }
  k = max(foldid)
  fold_sizes = table(foldid)
  fold_weights = fold_sizes/sum(fold_sizes)

  for (i in 1:(nrow(cross_validation_grid))) {
    fold_performance_measure = array(NA, dim = c(k, length(nprune_values)))
    num_intervals = cross_validation_grid$num_bins[i]
    interval = seq(0, time_point, length.out = num_intervals + 1)
    degree = cross_validation_grid$degree_values[i]
    for (j in 1:k) {
      train_data = data[foldid != j, ]
      test_data = data[foldid == j, ]
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
      # Get the expanded train dataset using combined IPCW
      expanded_train_data_with_weight = bin_combined_ipcw(interval = interval,data = train_data,
                                                          var_name = var_name, time_point = time_point,
                                                          return_type="dataset")
      valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
      filtered_train_E = expanded_train_data_with_weight$E[valid_indices]
      filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
      # Fit MARS with BINNED IPCW weights
      # rss_vals = numeric(length(var_name))
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
                          nprune = NULL,
                          degree = degree,
                          ...)
      for(np in seq_along(nprune_values)){
        nested_nprune = nprune_values[np]
        if(is.na(nested_nprune)){
          nested_nprune = NULL
        }
        earth_prune_model = update(earth_model, nprune = nested_nprune)
        est_test_event_single = predict(earth_prune_model, newdata = test_data, type = "response")
        est_test_surv_probability = (1-est_test_event_single)^num_intervals
        est_test_event_probability = 1-est_test_surv_probability
        fold_performance_measure[j,np] <- switch(measure,
                                              "C_index" = concordance.index(x = est_test_event_probability,
                                                                            surv.time = test_data$observed_time,
                                                                            surv.event = test_data$sigma)$c.index,
                                              "Log_Likelihood_Neg" = weighted_loglikelihood(IPCW = test_IPCW,
                                                                                            status = test_data$E,
                                                                                            est_test_event_probability),
                                              "Brier_Score" = Weighted_Brier_Score(est_test_event_probability,
                                                                                   status = test_data$E,
                                                                                   IPCW = test_IPCW),
                                              stop("Invalid measure. Please use 'C_index', 'Log_Likelihood_Neg', or 'Brier_Score'.")
        )
      }
    }
    for(np in seq_along(nprune_values)){
      valid_indices <- !is.na(fold_performance_measure[, np])
      fold_weights = fold_sizes[valid_indices]/sum(fold_sizes[valid_indices])
      cv_results[[measure]][cv_results$num_bins == num_intervals & 
                              cv_results$degree_values == degree & 
                              cv_results$nprune_values == nprune_values[np]] = sum(fold_performance_measure[valid_indices, np] * fold_weights)
    }
  }
  return(list(results = cv_results))
}

