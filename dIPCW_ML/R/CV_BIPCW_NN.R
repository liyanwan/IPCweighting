## ---------------------------------------------------------------------------------------------------------------
#' Compute Cross-Validation Results Using an IPC weighted Neural Network Model with discrete-time modelling
#'
#' @param X A dataframe of covariates.
#' @param Y A dataframe containing observed_time, sigma, and E.
#' @param time_point The specified cut-off time for evaluation.
#' @param measure The metric used for selecting the optimal tuning parameters. 
#' @param alpha A numeric vector specifying momentum rate for GD to evaluate.
#' @param lambda A numeric vector specifying L2 regularization parameter to evaluate.
#' @param lr_rate A numeric vector specifying learning rate parameter to evaluate.
#' @param range_intervals A numeric vector specifying the number of intervals to evaluate during cross-validation.
#' @param k The number of cross-validation folds.
#' @param foldid A vector indicating the fold assignment for each observation.
#' @param proxy_data A dataframe that may be used to compute IPCW in the test set.
#' @param sample_batch If TRUE, ensures that the same set of samples is selected in each batch across cross-validation folds.
#' @param batch_seed Used with set.seed(batch_seed) to guarantee reproducibility of batch sampling when sample_batch = TRUE.
#' 
fit_deepnn_binnedweights <- function(X,
                                     Y,
                                     time_point,
                                     measure = c('C_index', 'Log_Likelihood_Neg', 'Brier_Score'),
                                     alpha,
                                     lambda,
                                     lr_rate,
                                     range_intervals, 
                                     k = 5, 
                                     foldid = NULL,
                                     proxy_data = NULL, 
                                     model,
                                     epochs = 200,
                                     batch_size = 64,
                                     epsilon = 1e-3,
                                     verbose = FALSE,
                                     sample_batch = TRUE,
                                     batch_seed = 0,
                                      ...) {
  data = cbind(X, Y)
  var_name = colnames(X)
  
  cv_results = subset(expand.grid(
    num_bins = range_intervals,
    alpha = alpha,
    lambda = lambda,
    lr_rate = lr_rate,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ), lr_rate*lambda <= 0.1)
  cv_results[[measure]] = NA
  
  if (is.null(foldid)) foldid = sample(rep(seq_len(k), length.out = nrow(X)))
  k = max(foldid)
  fold_sizes = as.integer(table(factor(foldid, levels = seq_len(k))))
  fold_weights = fold_sizes / sum(fold_sizes)
  
  for (i in seq_len(nrow(cv_results))) {
    print(i)
    num_intervals = cv_results$num_bins[i]
    interval = seq(0, time_point, length.out = num_intervals + 1)
    al = cv_results$alpha[i]
    lam = cv_results$lambda[i]
    lr = cv_results$lr_rate[i]
    fold_performance_measure = rep(NA, k)
    
    for (f in seq_len(k)) {
      idx_tr = which(foldid != f)
      idx_te = which(foldid == f)
      train_data = data[idx_tr, , drop = FALSE]
      test_data  = data[idx_te, , drop = FALSE]
      train_Ghat = Get_Ghat(observed_time = train_data[, 'observed_time'], sigma = train_data[,'sigma'])
      # Compute IPC weights of test data using G from train data
      if (!is.null(proxy_data)) {
        test_proxy = rbind(test_data, proxy_data)
        totaltest_G = Ghat_newtime(step_function = train_Ghat, 
                                   new_observed_time = test_proxy[, 'observed_time'], 
                                   time_point = time_point)
        total_IPCW = ifelse(is.na(test_data[,'E']), 0, 1/totaltest_G)
        test_IPCW = total_IPCW[1:nrow(test_data)]
      } else {
        test_G = Ghat_newtime(step_function = train_Ghat,
                              new_observed_time = test_data[, 'observed_time'],
                              time_point = time_point)
        test_IPCW = ifelse(is.na(test_data[,'E']), 0, 1/test_G)
      }
      expanded = bin_combined_ipcw(interval = interval, data = train_data, var_name = var_name,
                                    time_point = time_point, return_type = "dataset")
      valid_idx = which(is.finite(expanded$IPCW) & expanded$IPCW > 0)
      if (!length(valid_idx)) next
      model_data = expanded[valid_idx, , drop = FALSE]
      filtered_train_E = model_data$E
      filtered_weights = model_data$IPCW
      if(sample_batch) set.seed(batch_seed)
      fit = dnn::deepGlm(filtered_train_E ~ as.matrix(model_data[var_name]),
                          model = model, family = "binomial", weights = filtered_weights,
                          epochs = epochs, lr_rate = lr, batch_size = batch_size, 
                          alpha = al, lambda = lam, verbose = verbose)
      est_test_event_single = as.numeric(predict(fit$model, x = as.matrix(test_data[,var_name, drop=FALSE])))
      est_test_surv_probability = (1 - est_test_event_single)^num_intervals
      est_test_event_probability = 1 - est_test_surv_probability
      val = switch(
        measure,
        C_index = concordance.index(x = est_test_event_probability,
                                    surv.time = test_data[, 'observed_time'],
                                    surv.event = test_data[,'sigma'])$c.index,
        Log_Likelihood_Neg = weighted_loglikelihood(IPCW = test_IPCW,
                                                    status = test_data[,'E'],
                                                    event_prob = est_test_event_probability),
        Brier_Score = Weighted_Brier_Score(event_prob = est_test_event_probability,
                                           status = test_data[,'E'],
                                           IPCW = test_IPCW),
        stop("Invalid measure. Use 'C_index', 'Log_Likelihood_Neg', or 'Brier_Score'.")
      )
      if (is.finite(val)) fold_performance_measure[f] = val
    }
    ok = is.finite(fold_performance_measure)
    if (all(ok)) {
      cv_results[[measure]][i] = sum(fold_performance_measure[ok] * w)
    }
  }
  list(results = cv_results, measure = measure)
}