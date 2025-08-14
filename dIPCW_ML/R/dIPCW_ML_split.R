dIPCW_ML_split <- function(time_point,
                           X,
                           Y,
                           train_prop,
                           measure,
                           range_intervals,
                           true_surv = NULL,
                           learner_list,
                           params_list,
                           proxy_data = NULL,
                           learner_k = 5,
                           learner_foldid = NULL,
                           include_naive = FALSE,
                           include_SA = FALSE,
                           surv_params = list(),
                           include_opt = FALSE,
                           max_intervals = max(range_intervals),
                           min_intervals = min(range_intervals)) {
  X = data.frame(X)
  Y = data.frame(Y)
  E_factor = factor(Y$E, levels = c(0, 1, NA), exclude = NULL)
  train_index <- createDataPartition(E_factor, p = train_prop, list = FALSE)
  M = Y$M
  event_time = Y$event_time
  if("M" %in% colnames(Y)){
    Y = within(Y, rm("M"))
  }
  if("event_time" %in% colnames(Y)){
    Y = within(Y, rm("event_time"))
  }
  train_X <- X[train_index, ]
  test_X <- X[-train_index, ]
  train_Y <- Y[train_index, ]
  test_Y <- Y[-train_index, ]
  test_M <- M[-train_index]
  test_event_time = event_time[-train_index]
  train_data = cbind.data.frame(train_X, train_Y)
  test_data = cbind.data.frame(test_X, test_Y)
  var_name = colnames(X)
  if(!is.null(true_surv)){
    test_true_surv = true_surv[-train_index]
  }
  if(is.null(learner_foldid)){
    trainE_stratified = ifelse(is.na(train_Y$E), "missing", as.character(train_Y$E))
    train_E_factor = factor(trainE_stratified)
    cv_folds <- createFolds(train_E_factor, k = learner_k, returnTrain = FALSE)
    learner_foldid = rep(NA, nrow(train_Y))
    for (i in seq_along(cv_folds)) {
      learner_foldid[cv_folds[[i]]] = i
    }
  }
  train_Ghat = Get_Ghat(observed_time = train_data$observed_time, sigma = train_data$sigma)
  test_G = Ghat_newtime(step_function = train_Ghat,
                        new_observed_time = test_data$observed_time,
                        time_point = time_point)
  test_IPCW = ifelse(is.na(test_data$E), 0, 1/test_G)
  # naive learners
  if(include_naive){
    naive_indices = which(!is.na(train_data$E))
    naive_train_data = train_data[naive_indices, ]
    naive_train_X = naive_train_data[var_name]
    naive_train_Y = data.frame(E = naive_train_data$E, sigma = naive_train_data$sigma, observed_time = naive_train_data$observed_time)
    naive_learner_foldid =learner_foldid[naive_indices]
    naive_learners_test_EP = fit_singlebin(train_X = naive_train_X, train_Y = naive_train_Y, X_validation = test_X, time_point = time_point,
                                           measure = measure, learner_list, params_list, return_type = "naive", 
                                           k = learner_k, foldid = naive_learner_foldid, proxy_data = proxy_data)$single_test_event_probability
  }
  # Extract test EP from Binned version learners
  model_list_with_bin <- train_all_learners_all(train_X = train_X, train_Y = train_Y, X_validation = test_X, 
                                                time_point = time_point, measure = measure, learner_list = learner_list, 
                                                params_list = params_list, range_intervals = range_intervals, k = learner_k, 
                                                proxy_data = proxy_data, foldid = learner_foldid)
  if(include_opt){
    model_list_opt <- train_all_learners_all(train_X = train_X, train_Y = train_Y, X_validation = test_X, 
                                             time_point = time_point, measure = measure, learner_list, params_list,
                                             range_intervals = c(min_intervals: max_intervals), k = learner_k, 
                                             proxy_data = proxy_data, foldid = learner_foldid)
  } else {
    model_list_opt = NULL
  }
  
  # Coxph Survival Model
  if(include_SA){
    if(any(sapply(learner_list, identical, y = base_mars_BinnedIPCW))){
      add_cox_test_EP = additive_cox(train_data = train_data, test_X = test_X, time_point = time_point,
                                     cts.num = ifelse(!is.null(surv_params$cts.num), surv_params$cts.num, 5),
                                     k = ifelse(!is.null(surv_params$k), surv_params$k, 10))
    }
    if(any(sapply(learner_list, identical, y = base_glmnet_BinnedIPCW))){
      index_match <- which(sapply(learner_list, identical, y = base_glmnet_BinnedIPCW))
      cox_penalized_test_EP = do.call(cbind,
                                      lapply(index_match, function(index_cox){
                                        coxph_with_penalty(
                                          train_data = train_data,
                                          test_X = test_X,
                                          time_point = time_point,
                                          alpha_num = params_list[[index_cox]]$alpha_num,
                                          measure = measure,
                                          lambda_grid = params_list[[index_cox]]$lambda_values,
                                          useMin = params_list[[index_cox]]$useMin,
                                          foldid = learner_foldid
                                        )
                                      })
      )
    }
    if(any(sapply(learner_list, identical, y = base_glm_BinnedIPCW))){
      index_match <- which(sapply(learner_list, identical, y = base_glm_BinnedIPCW))
      cox_penalized_test_EP = do.call(cbind,
                                      lapply(index_match, function(index_cox){
                                        coxph_with_penalty(
                                          train_data = train_data,
                                          test_X = test_X,
                                          time_point = time_point,
                                          alpha_num = NULL,
                                          measure = measure,
                                          lambda_grid = NULL,
                                          useMin = params_list[[index_cox]]$useMin
                                        )
                                      })
      )
    }
    # Survival Tree model (Assume Poisson Distribution (exp form))
    if(any(sapply(learner_list, identical, y = base_CTree_BinnedIPCW))){
      formula = as.formula(paste("Surv(observed_time, sigma) ~", paste(var_name, collapse = " + ")))
      surv_tree = rpart(formula, data = train_data)
      best_cp <- surv_tree$cptable[which.min(surv_tree$cptable[, "xerror"]), "CP"]
      pruned_tree <- prune(surv_tree, cp = best_cp)
      hazards = predict(pruned_tree, newdata = test_data, type = "matrix")[, 1]
      s0 = survreg(Surv(observed_time, sigma) ~ 1, data = train_data, dist = "exponential")
      e0 = exp(-summary(s0)$coefficients[1])
      surv_tree_surv_prob = exp(-e0 * hazards * time_point)
      surv_tree_event_probability = 1 - surv_tree_surv_prob
      STree_test_EP = data.frame(SurvTree = surv_tree_event_probability)
    }
  }
  all_model_event_probability_layer <- data.frame()
  all_model_optimal_bin = c()
  for(i in 1:length(model_list_with_bin)){
    item_bin = model_list_with_bin[[i]]$bin
    class_name = model_list_with_bin[[i]]$class_name
    pred_test_event_probability = as.vector(model_list_with_bin[[i]]$test_event_pred)
    col_name = paste0(class_name, item_bin)
    if(nrow(all_model_event_probability_layer)==0){
      all_model_event_probability_layer = data.frame(pred_test_event_probability)
      all_model_optimal_bin = c(item_bin)
      names(all_model_event_probability_layer) = col_name
    } else{
      all_model_event_probability_layer[[col_name]] = pred_test_event_probability
      all_model_optimal_bin = c(all_model_optimal_bin, item_bin)
    }
  }
  if(include_opt){
    for(i in 1:length(model_list_opt)){
      item_bin = model_list_opt[[i]]$bin
      class_name = model_list_opt[[i]]$class_name
      pred_test_event_probability = as.vector(model_list_opt[[i]]$test_event_pred)
      col_name = paste0("OPT", class_name)
      if(nrow(all_model_event_probability_layer)==0){
        all_model_event_probability_layer = data.frame(pred_test_event_probability)
        all_model_optimal_bin = c(item_bin)
        names(all_model_event_probability_layer) = col_name
      }
      else{
        all_model_event_probability_layer[[col_name]] = pred_test_event_probability
        all_model_optimal_bin = c(all_model_optimal_bin, item_bin)
      }
    }
  }
  if(include_naive){
    all_model_event_probability_layer = cbind(all_model_event_probability_layer,
                                              naive_learners_test_EP)
  }
  if(include_SA){
    if (exists("cox_penalized_test_EP")) {
      all_model_event_probability_layer <- cbind(all_model_event_probability_layer,
                                                 cox_penalized_test_EP)
    }
    if (exists("cox_test_EP")) {
      all_model_event_probability_layer <- cbind(all_model_event_probability_layer,
                                                 cox_test_EP)
    }
    if(exists("add_cox_test_EP")){
      all_model_event_probability_layer <- cbind(all_model_event_probability_layer,
                                                 add_cox_test_EP)
    }
    if (exists("STree_test_EP")) {
      all_model_event_probability_layer = cbind(all_model_event_probability_layer, STree_test_EP)
    }
  }
  measure_results = vector("list", ncol(all_model_event_probability_layer))
  ols_layer = vector("list", ncol(all_model_event_probability_layer))
  C_layer = vector("list", ncol(all_model_event_probability_layer))
  C_true_layer = vector("list", ncol(all_model_event_probability_layer))
  AUC_layer = vector("list", ncol(all_model_event_probability_layer))
  NCLL_layer = vector("list", ncol(all_model_event_probability_layer))
  NCBS_layer = vector("list", ncol(all_model_event_probability_layer))
  for (kk in seq_len(ncol(all_model_event_probability_layer))) {
    predictions <- all_model_event_probability_layer[, kk]
    C_layer[[kk]] = concordance.index(x = predictions,
                                      surv.time = test_Y$observed_time,
                                      surv.event = test_Y$sigma)$c.index
    C_true_layer[[kk]] = concordance.index(x = predictions, 
                                           surv.time = test_event_time, 
                                           surv.event = rep(1, length(test_Y$sigma)))$c.index
    AUC_layer[[kk]] = auc(roc(test_M, predictions, quiet = TRUE))
    NCLL_layer[[kk]] = weighted_loglikelihood(rep(1, length(test_M)), test_M, predictions)
    NCBS_layer[[kk]] = Weighted_Brier_Score(predictions, test_M, rep(1, length(test_M)))
    if(!is.null(true_surv)){
      ols_layer[[kk]] = ols_error(predictions, 1 - test_true_surv)
    }
  }
  if(all(sapply(ols_layer, is.null))){
    results_df <- setNames(
      data.frame(
        Method = colnames(all_model_event_probability_layer),
        C = unlist(C_layer),
        C_2 = unlist(C_true_layer),
        AUC = unlist(AUC_layer),
        LL = unlist(NCLL_layer),
        BS = unlist(NCBS_layer),
        Opt_bins = c(all_model_optimal_bin, rep(NA, ncol(all_model_event_probability_layer) - length(all_model_optimal_bin)))
      ), c("Method","C-observed", "C-event", "AUC", "-LogL", "BS", "Opt_bins"))
  }else{
    results_df <- setNames(
      data.frame(
        Method = colnames(all_model_event_probability_layer),
        C = unlist(C_layer),
        C_2 = unlist(C_true_layer),
        AUC = unlist(AUC_layer),
        LL = unlist(NCLL_layer),
        BS = unlist(NCBS_layer),
        OLS = unlist(ols_layer),
        Opt_bins = c(all_model_optimal_bin, rep(NA, ncol(all_model_event_probability_layer) - length(all_model_optimal_bin)))
      ), c("Method","C-observed", "C-event", "AUC", "-LogL", "BS", "OLS", "Opt_bins"))
  }
  return(list(results_df = results_df,
              all_layer_testEP = all_model_event_probability_layer,
              model_list_with_bin = model_list_with_bin,
              model_list_opt = model_list_opt))
}

