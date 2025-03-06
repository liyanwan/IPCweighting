dIPCW_SL<- function(time_point,
                   X,
                   Y,
                    newX = NULL,
                    newY = NULL,
                   measure,
                   range_intervals,
                   true_surv = NULL,
                   SL_library,
                   meta_learner_params = NULL,
                   method,
                   control = list(),
                   cvControl = list(),
                   verbose = FALSE) {
  if(is.null(newX)){
    newX = X
  }
  if(is.null(newY)){
    newY = Y
  }
  X = data.frame(X)
  Y = data.frame(Y)
  newX = data.frame(newX)
  newY = data.frame(newY)
  train_data = cbind.data.frame(X, Y)
  test_data = cbind.data.frame(newX, newY)
  var_name = colnames(X)
  if(!is.null(true_surv)){
    test_true_surv = true_surv[-train_index]
  }
  train_Ghat = Get_Ghat(observed_time = train_data$observed_time, sigma = train_data$sigma)
  test_G = Ghat_newtime(step_function = train_Ghat,
                        new_observed_time = test_data$observed_time,
                        time_point = time_point)
  test_IPCW = ifelse(test_data$sigma == 0 & test_data$observed_time<time_point, 0, 1/test_G)
  SL_binned_info = dIPCWSuperLearner(Y = Y,
                                     X = X,
                                     time_point = time_point,
                                     newX = newX,
                                     family = binomial(),
                                     SL.library = SL_library,
                                     measure = measure,
                                     method = method,
                                     id = NULL,
                                     verbose = verbose,
                                     control = control,
                                     meta_learner_params = meta_learner_params,
                                     cvControl = cvControl,
                                     env = parent.frame())
  weighted_event_probability_SL_binned = SL_binned_info$SL.predict
  coef = SL_binned_info$coef
  ML_CVRisk = SL_binned_info$MLCVRisk
  if(!is.null(ML_CVRisk)){
    names(ML_CVRisk) = measure
  }
  head_name = sub("method\\.", "", method)
  if((length(range_intervals) == 1)&(!is.null(meta_learner_params))){
    SL_colnames = paste0(head_name, "SLIPCW_", range_intervals)
  } else {
    SL_colnames = paste0(head_name, "BinnedSL_IPCW")
  }
  all_model_event_probability_layer <- data.frame()
  all_model_optimal_bin = c()
  colnames(weighted_event_probability_SL_binned) = SL_colnames
  all_model_event_probability_layer = data.frame(weighted_event_probability_SL_binned)
  measure_results = vector("list", ncol(all_model_event_probability_layer))
  ols_layer = vector("list", ncol(all_model_event_probability_layer))
  C_layer = vector("list", ncol(all_model_event_probability_layer))
  LL_layer = vector("list", ncol(all_model_event_probability_layer))
  BS_layer = vector("list", ncol(all_model_event_probability_layer))
  AUC_layer = vector("list", ncol(all_model_event_probability_layer))
  
  for (kk in seq_len(ncol(all_model_event_probability_layer))) {
    predictions <- all_model_event_probability_layer[, kk]
    C_layer[[kk]] = concordance.index(x = predictions,
                                      surv.time = newY$observed_time,
                                      surv.event = newY$sigma)$c.index
    LL_layer[[kk]] = weighted_loglikelihood(test_IPCW, newY$E, predictions)
    BS_layer[[kk]] = Weighted_Brier_Score(predictions, newY$E, test_IPCW)
    non_censored_ind = !is.na(newY$E)
    AUC_layer[[kk]] = auc(roc(newY$E[non_censored_ind], predictions[non_censored_ind]))
    
    if(!is.null(true_surv)){
      ols_layer[[kk]] = ols_error(predictions, 1 - test_true_surv)
    }
  }
  if(all(sapply(ols_layer, is.null))){
    results_df <- setNames(
      data.frame(
        Method = colnames(all_model_event_probability_layer),
        C = unlist(C_layer),
        LL = unlist(LL_layer),
        BS = unlist(BS_layer),
        AUC = unlist(AUC_layer),
        Opt_bins = rep(NA, ncol(all_model_event_probability_layer))
      ), c("Method","C-index", "-Log-Likelihood", "Brier_Score", "AUC", "Opt_bins"))
  }else{
    results_df <- setNames(
      data.frame(
        Method = colnames(all_model_event_probability_layer),
        C = unlist(C_layer),
        LL = unlist(LL_layer),
        BS = unlist(BS_layer),
        AUC = unlist(AUC_layer),
        OLS = unlist(ols_layer),
        Opt_bins = rep(NA, ncol(all_model_event_probability_layer))
      ), c("Method","C-index", "-Log-Likelihood", "Brier_Score", "AUC", "OLS_error", "Opt_bins"))
  }
  return(list(results_df = results_df,
              all_layer_testEP = all_model_event_probability_layer,
              coef = coef,
              ML_CVRisk=ML_CVRisk))
}