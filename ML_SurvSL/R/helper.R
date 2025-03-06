## ---------------------------------------------------------------------------------------------------------------
Get_Ghat <-function(observed_time, sigma){
  # A helper function, using Kaplan Meier estimator of survival distribution of the censoring times.
  full_dt = data.frame(observed_time, sigma)
  mm = subjectWeights(formula = Surv(observed_time, sigma)~1,
                      data = full_dt,
                      method = "marginal")
  step_function <- function(time) {
    predict(mm$fit, newdata = full_dt, 
            times = time - min(full_dt$observed_time) * 1e-5, 
            level.chaos = 1, mode = "matrix", type = "surv")
  }
  return(step_function)
}


## ---------------------------------------------------------------------------------------------------------------
Ghat_newtime <- function(step_function, new_observed_time, time_point){
  ordered_weights = step_function(new_observed_time)
  new_G = ordered_weights[rank(new_observed_time)]
  new_G[new_observed_time >= time_point] = step_function(time_point)
  return(new_G)
}

## ---------------------------------------------------------------------------------------------------------------
get_IPCW_relevant <- function(observed_time, sigma, time_point, return_type = "data"){
  # km_censor_Xi = survfit(Surv(observed_time, 1-sigma)~1, data = dt)
  # survest_Xi = stepfun(km_censor_Xi$time, c(1, km_censor_Xi$surv))
  # censor_prob_Xi = survest_Xi(ifelse(dt$observed_time<time_point, dt$observed_time,time_point))
  # dt$G_hat_Vi = censor_prob_Xi
  dt = data.frame(observed_time, sigma)
  mm = subjectWeights(formula = Surv(observed_time, sigma)~1,
                      data = dt,
                      method = "marginal")
  time_pointG = predict(mm$fit,newdata=dt,times=time_point,level.chaos=1,mode="matrix",type="surv")
  dt$G_hat_Vi <- mm$weights[rank(dt$observed_time)]
  dt$G_hat_Vi[dt$observed_time >= time_point] = time_pointG
  IPCW = ifelse(dt$sigma==0 & dt$observed_time<=time_point, 0, 1/pmax(dt$G_hat_Vi, 1e-15))
  if (return_type == "IPCW"){
    return(IPCW)
  }
  else{
    dt$IPCW = IPCW
    return(dt)
  }
}


## ---------------------------------------------------------------------------------------------------------------
bin_combined_ipcw_separate <- function(interval, X, Y, time_point, return_type="model"){
  total=data.frame()
  if (is.matrix(X)) {
    X <- as.data.frame(X)
  }
  if (is.matrix(Y)) {
    Y <- as.data.frame(Y)
  }
  train_data = cbind(X, Y)
  var_name = colnames(X)
  Ghat = Get_Ghat(observed_time = Y$observed_time, sigma = Y$sigma)
  for(i in 1:(length(interval)-1)){
    start = interval[i]
    end = interval[i+1]
    numerator = Ghat(start)
    batch = train_data[train_data$observed_time>start, ]
    batch$E = get_status(batch$observed_time, batch$sigma, end)
    batch$G = Ghat_newtime(step_function = Ghat, new_observed_time = batch$observed_time, time_point = end)
    batch$IPCW = ifelse(batch$sigma == 0 & batch$observed_time<=end, 0, numerator/batch$G)
    if (is.null(total)) {
      total = data.frame(batch)
    } else {
      total = rbind(total, batch)
    }
  }
  if(return_type == "dataset"){
    return(total)
  }
  else if(return_type == "model"){
    formula = as.formula(paste("E ~", paste(var_name, collapse = " + ")))
    model = glm(formula, data = total, family = binomial, weights = IPCW)
    return(model)
  }
}


## ---------------------------------------------------------------------------------------------------------------
cross_validate_lambda_SingleMeasure <- function(X, Y, time_point, alpha_num, lambda_values, measure,
                                                range_intervals = c(1:60), k = 5, proxy_data = NULL, foldid = NULL, ...) {
  var_name = colnames(X)
  n = nrow(X)
  nlambda = length(lambda_values)
  cv_results = data.frame(num_bins = range_intervals,
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
      train_X = X[foldid != f, ]
      test_X = X[foldid == f, ]
      train_Y = Y[foldid != f, ]
      test_Y = Y[foldid == f, ]
      train_Ghat = Get_Ghat(observed_time = train_Y$observed_time, sigma = train_Y$sigma)
      # Compute IPC weights of test data using G from train data
      if (!is.null(proxy_data)) {
        test_proxy = rbind(test_Y, proxy_data)
        totaltest_G = Ghat_newtime(step_function = train_Ghat, 
                                   new_observed_time = test_proxy$observed_time, 
                                   time_point = time_point)
        total_IPCW = ifelse(test_Y$sigma == 0 & test_Y$observed_time<=time_point, 0, 1/totaltest_G)
        test_IPCW = total_IPCW[1:nrow(test_Y)]
      } else {
        test_G = Ghat_newtime(step_function = train_Ghat,
                              new_observed_time = test_Y$observed_time,
                              time_point = time_point)
        test_IPCW = ifelse(test_Y$sigma == 0 & test_Y$observed_time<=time_point, 0, 1/test_G)
      }
      expanded_train_data_with_weight = bin_combined_ipcw_separate(interval = interval, 
                                                                   X = train_X,
                                                                   Y = train_Y,
                                                                   time_point = time_point, 
                                                                   return_type = "dataset")
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
          standardize = TRUE
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
                                                                                    surv.time = test_Y$observed_time,
                                                                                    surv.event = test_Y$sigma)$c.index,
                                                   "Log_Likelihood_Neg" = weighted_loglikelihood(IPCW = test_IPCW, 
                                                                                                 status = test_Y$E, 
                                                                                                 event_prob = est_test_event_probability[, j]),
                                                   "Brier_Score" = Weighted_Brier_Score(event_prob = est_test_event_probability[, j], 
                                                                                        status = test_Y$E, 
                                                                                        IPCW = test_IPCW),
                                                   stop("Invalid measure. Please use 'C_index', 'Log_Likelihood_Neg', or 'Brier_Score'.")
          )
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


## --------------------------------------------------------------------------------------------------------------
get_optimal_row_indices_SingleMeasure <- function(optimal_results, useMin = TRUE, measure) {
  ncolumn = ncol(optimal_results)
  if(ncol(optimal_results) == 2){
    if(useMin == TRUE){
      if(measure == "C_index"){
        optimal_value = max(optimal_results[, ncolumn], na.rm = TRUE)
      } else {
        optimal_value = min(optimal_results[, ncolumn], na.rm = TRUE)
      }
      optimal_indices = which(optimal_results[, ncolumn] == optimal_value)
    }
    else{
      se <- sd(optimal_results[, ncolumn], na.rm = TRUE) / sqrt(sum(!is.na(optimal_results[, ncolumn])))
      if(measure == "C_index"){
        optimal_value = max(optimal_results[, ncolumn], na.rm = TRUE)
        optimal_indices = which(optimal_results[, ncolumn] >= (optimal_value - se))
      } else {
        optimal_value = min(optimal_results[, ncolumn], na.rm = TRUE)
        optimal_indices = which(optimal_results[, ncolumn] <= (optimal_value + se))
      }
    }
    return(list(optimal_index = min(optimal_indices)))
  }
  else if (ncolumn == 3){
    if(useMin == TRUE){
      if(measure == "C_index"){
        optimal_value = max(optimal_results[, ncolumn], na.rm = TRUE)
      } else {
        optimal_value = min(optimal_results[, ncolumn], na.rm = TRUE)
      }
      optimal_indices = which(optimal_results[, ncolumn] == optimal_value)
      get_best_index <- function(optimal_results, indices) {
        candidates <- optimal_results[indices, ]
        if (nrow(candidates) == 0) {
          return(NULL)
        }
        min_value_col1 <- min(optimal_results[, 1], na.rm = TRUE)
        max_value_col2 <- max(optimal_results[, 2], na.rm = TRUE)
        candidates <- candidates %>%
          mutate(score = (min_value_col1 / .[, 1]) + (.[, 2] / max_value_col2)) %>%
          arrange(desc(score)) %>%
          head(1)
        return(which(optimal_results[, 1] == candidates[, 1] & optimal_results[, 2] == candidates[, 2]))
      }
      best_index <- get_best_index(optimal_results, optimal_indices)
    }
    else{
      se <- sd(optimal_results[, ncolumn], na.rm = TRUE) / sqrt(sum(!is.na(optimal_results[, ncolumn])))
      if(measure == "C_index"){
        optimal_value = max(optimal_results[, ncolumn], na.rm = TRUE)
        optimal_indices = which(optimal_results[, ncolumn] >= (optimal_value - se))
      } else {
        optimal_value = min(optimal_results[, ncolumn], na.rm = TRUE)
        optimal_indices = which(optimal_results[, ncolumn] <= (optimal_value + se))
      }
      get_best_index <- function(optimal_results, indices) {
        candidates <- optimal_results[indices, ]
        if (nrow(candidates) == 0) {
          return(NULL)
        }
        min_value_col1 <- min(optimal_results[, 1], na.rm = TRUE)
        max_value_col2 <- max(optimal_results[, 2], na.rm = TRUE)
        candidates <- candidates %>%
          mutate(score = (min_value_col1 / .[, 1]) + (.[, 2] / max_value_col2)) %>%
          arrange(desc(score)) %>%
          head(1)
        return(which(optimal_results[, 1] == candidates[, 1] & optimal_results[, 2] == candidates[, 2]))
      }
      best_index <- get_best_index(optimal_results, optimal_indices)
    }
    return(list(optimal_index = best_index))
  }
  else{
    if(useMin == TRUE){
      if(measure == "C_index"){
        optimal_value = max(optimal_results[, ncolumn], na.rm = TRUE)
      } else {
        optimal_value = min(optimal_results[, ncolumn], na.rm = TRUE)
      }
      optimal_indices = which(optimal_results[, ncolumn] == optimal_value)
      get_best_index <- function(optimal_results, indices) {
        candidates <- optimal_results[indices, ]
        if (nrow(candidates) == 0) {
          return(NULL)
        }
        min_value_col1 <- min(optimal_results[, 1], na.rm = TRUE)
        min_value_col2 <- min(optimal_results[, 2], na.rm = TRUE)
        min_value_col3 <- min(optimal_results[, 3], na.rm = TRUE)
        candidates <- candidates %>%
          mutate(score = (min_value_col1 / .[, 1]) + (min_value_col2 / .[, 2]) + (min_value_col3 / .[, 3])) %>%
          arrange(desc(score)) %>%
          head(1)
        return(which(optimal_results[, 1] == candidates[, 1] & optimal_results[, 2] == candidates[, 2] & optimal_results[, 3] == candidates[, 3]))
      }
      best_index <- get_best_index(optimal_results, optimal_indices)
    }
    else{
      se <- sd(optimal_results[, ncolumn], na.rm = TRUE) / sqrt(sum(!is.na(optimal_results[, ncolumn])))
      if(measure == "C_index"){
        optimal_value = max(optimal_results[, ncolumn], na.rm = TRUE)
        optimal_indices = which(optimal_results[, ncolumn] >= (optimal_value - se))
      } else {
        optimal_value = min(optimal_results[, ncolumn], na.rm = TRUE)
        optimal_indices = which(optimal_results[, ncolumn] <= (optimal_value + se))
      }
      get_best_index <- function(optimal_results, indices) {
        candidates <- optimal_results[indices, ]
        if (nrow(candidates) == 0) {
          return(NULL)
        }
        min_value_col1 <- min(optimal_results[, 1], na.rm = TRUE)
        min_value_col2 <- min(optimal_results[, 2], na.rm = TRUE)
        min_value_col3 <- min(optimal_results[, 3], na.rm = TRUE)
        candidates <- candidates %>%
          mutate(score = (min_value_col1 / .[, 1]) + (min_value_col2 / .[, 2]) + (min_value_col3 / .[, 3])) %>%
          arrange(desc(score)) %>%
          head(1)
        return(which(optimal_results[, 1] == candidates[, 1] & optimal_results[, 2] == candidates[, 2] & optimal_results[, 3] == candidates[, 3]))
      }
      best_index <- get_best_index(optimal_results, optimal_indices)
    }
    return(list(optimal_index = best_index))
  }
}
