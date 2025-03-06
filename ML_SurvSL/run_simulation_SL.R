run_simulation_SL <- function(num_obs,
                              time_point,
                              X,
                              lc,
                              dist,
                              params,
                              measure,
                              censor_dist,
                              censor_params,
                              range_intervals,
                              train_prop,
                              true_surv = NULL,
                              SL_library,
                              meta_learner_params = NULL,
                              method,
                              control = list(),
                              cvControl = list(),
                              verbose = FALSE) {
  params$lc <- lc
  power <- exp(params$lc)
  dt <- simulation_data(num_obs, dist, params, censor_dist, censor_params, time_point, X)
  dt$power <- power
  X <- as.data.frame(X)
  train_index <- sample(seq_len(nrow(dt)), size = floor(train_prop * nrow(dt)))
  train_data <- dt[train_index, ]
  test_data <- dt[-train_index, ]
  train_X <- X[train_index, ]
  test_X <- X[-train_index, ]
  # Y include E, censored status and observed_time
  Y <- data.frame(E = dt$E, sigma = dt$sigma, observed_time = dt$observed_time)
  train_Y <- Y[train_index, ]
  test_Y <- Y[-train_index, ]
  if(!is.null(true_surv)){
    test_true_surv = true_surv[-train_index]
  }
  var_name = colnames(X)
  train_Ghat = Get_Ghat(observed_time = train_data$observed_time, sigma = train_data$sigma)
  test_G = Ghat_newtime(step_function = train_Ghat,
                        new_observed_time = test_data$observed_time,
                        time_point = time_point)
  test_IPCW = ifelse(test_data$sigma == 0 & test_data$observed_time<time_point, 0, 1/test_G)
  SL_binned_info = dIPCWSuperLearner(Y = train_Y,
                                     X = train_X,
                                     time_point = time_point,
                                     newX = test_X,
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
  for (kk in seq_len(ncol(all_model_event_probability_layer))) {
    predictions <- all_model_event_probability_layer[, kk]
    if (measure == "C_index") {
      measure_results[[kk]] <- concordance.index(x = predictions,
                                                 surv.time = test_Y$observed_time,
                                                 surv.event = test_Y$sigma)$c.index
    } else if (measure == "Log_Likelihood_Neg") {
      measure_results[[kk]] <- weighted_loglikelihood(test_IPCW, test_Y$E, predictions)
    } else if (measure == "Brier_Score") {
      measure_results[[kk]] <- Weighted_Brier_Score(predictions, test_Y$E, test_IPCW)
    }
    C_layer[[kk]] = concordance.index(x = predictions,
                                      surv.time = test_Y$observed_time,
                                      surv.event = test_Y$sigma)$c.index
    if(!is.null(true_surv)){
      ols_layer[[kk]] = ols_error(predictions, 1 - test_true_surv)
    }
  }
  if(all(sapply(ols_layer, is.null))){
    results_df <- setNames(
      data.frame(
        Method = colnames(all_model_event_probability_layer),
        Value = unlist(measure_results),
        C = unlist(C_layer),
        Opt_bins = c(all_model_optimal_bin, rep(NA, ncol(all_model_event_probability_layer) - length(all_model_optimal_bin)))
      ), c("Method", measure, "C-index", "Opt_bins"))
  }else{
    results_df <- setNames(
      data.frame(
        Method = colnames(all_model_event_probability_layer),
        Value = unlist(measure_results),
        C = unlist(C_layer),
        OLS = unlist(ols_layer),
        Opt_bins = c(all_model_optimal_bin, rep(NA, ncol(all_model_event_probability_layer) - length(all_model_optimal_bin)))
      ), c("Method", measure, "C-index", "OLS_error", "Opt_bins"))
  }
  return(list(results_df = results_df,
              all_layer_testEP = all_model_event_probability_layer,
              coef = coef,
              ML_CVRisk=ML_CVRisk))
}
