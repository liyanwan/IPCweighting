## ---------------------------------------------------------------------------------------------------------------
#' Build a Suvival analysis model and estimate event probabilities on test_X.
#'
#' @param train_data A dataframe containing covariates observed_time, sigma, and E.
#' @param test_X A dataframe of covariates used for validation.
#' @param time_point The specified cut-off time for evaluation.
#' @param measure The metric used for selecting the optimal tuning parameters.
#' @param alpha_num If NULL, a standard Cox proportional hazards model (coxph) is used; 
#' otherwise, it represents the Elastic Net mixing parameter, ranging from 0 to 1.
#' @param lambda_grid If NULL, coxph is used; 
#' otherwise, it specifies a vector of lambda values for cross-validation in the Elastic Net model. 
#' @param useMin Logical; if TRUE, selects the tuning parameter based on the optimal measure value.
#' If FALSE, selects the parameter within one standard error of the optimal measure.

coxph_with_penalty <- function(train_data, test_X, time_point, measure = NULL, alpha_num = NULL, lambda_grid = NULL, foldid = NULL, useMin){
  var_name = colnames(test_X)
  formula = as.formula(paste("Surv(observed_time, sigma) ~", paste(var_name, collapse = " + ")))
  if(is.null(alpha_num)){
    cox_model = coxph(formula, data = train_data, method="breslow")
    lp = as.vector(rowSums(sweep(test_X, 2, as.vector(cox_model$coefficients),"*")))
    baseline_hazard = basehaz(cox_model,centered=FALSE)
    BLH_timepoint = baseline_hazard$hazard[which.min(abs(baseline_hazard$time - time_point))]
    coxph_survival_probability = exp(-BLH_timepoint*exp(lp))
    coxph_event_probability = 1 - coxph_survival_probability
    penalized_event_prob = data.frame(coxPH = coxph_event_probability)
    return(penalized_event_prob)
  }
  combined_penalized_event_prob = data.frame()
  Y_glmnet = Surv(train_data$observed_time, train_data$sigma)
  train_X = train_data[var_name]
  for (alpha in alpha_num){
    if (measure == "C_index"){
      model = cv.glmnet(as.matrix(train_X),
                        Y_glmnet,
                        family = "cox",
                        alpha = alpha,
                        lambda = lambda_grid,
                        type.measure = "C",
                        foldid = foldid)
    } else if (measure %in% c("Log_Likelihood_Neg", "Brier_Score")){
      model = cv.glmnet(as.matrix(train_X),
                        Y_glmnet,
                        family = "cox",
                        alpha = alpha,
                        lambda = lambda_grid,
                        type.measure = "deviance",
                        foldid = foldid)
    } else {
      stop("Invalid measure. Please use 'C_index', 'Log_Likelihood_Neg', or 'Brier_Score'.")
    }
    test_survfit <- survfit(formula = model, 
                            s = ifelse(useMin, "lambda.min", "lambda.1se"), 
                            x = as.matrix(train_X), 
                            y = Y_glmnet, 
                            newx = as.matrix(test_X))
    time_points = test_survfit$time
    index = max(which(time_points <= time_point))
    coxph_penal_survival_probability = test_survfit$surv[index, ]
    coxph_penal_event_probability = 1 - coxph_penal_survival_probability
    penalized_type = switch(
      as.character(alpha),
      "0" = "CoxPH Ridge",
      "1" = "CoxPH Lasso",
      paste0("CoxPH Elastic", alpha)
    )
    if(ncol(combined_penalized_event_prob) == 0){
      combined_penalized_event_prob = data.frame(coxph_penal_event_probability)
      colnames(combined_penalized_event_prob) = penalized_type
    }
    else{
      combined_penalized_event_prob[[penalized_type]] = coxph_penal_event_probability
    }
  }
  return(combined_penalized_event_prob)
}

## ---------------------------------------------------------------------------------------------------------------
#' Build an Additive Cox model and estimate event probabilities on test_X.
additive_cox <- function(train_data, test_X, time_point, cts.num=5, k = 10, var_threshold = ncol(test_X)){
  var_name = colnames(test_X)
  combined_gam_test_EP = list()
  train_X = train_data[var_name]
  cts.x <- apply(train_X, 2, function(x) (length(unique(x)) > cts.num))
  if (sum(!cts.x) > 0) {
    gam.model <- as.formula(paste("observed_time~", paste(paste("s(",
                                                                colnames(train_X[, cts.x, drop = FALSE]),
                                                                ", k = ", k, ")", sep = ""), collapse = "+"), "+", paste(
                                                                  colnames(train_X[,!cts.x, drop = FALSE]), collapse = "+")))
  } else {
    filtered_X = train_X[, cts.x, drop = FALSE]
    expanded_train_data_with_weight = bin_combined_ipcw(interval = c(0, time_point), data = train_data,
                                                        var_name = var_name, time_point = time_point, 
                                                        return_type="dataset")
    expanded_X = expanded_train_data_with_weight[var_name]
    valid_indices = which(expanded_train_data_with_weight$IPCW > 0)
    filtered_train_E = expanded_train_data_with_weight$E[valid_indices]
    filtered_train_X = expanded_X[valid_indices, ]
    filtered_weights = expanded_train_data_with_weight$IPCW[valid_indices]
    model_data = data.frame(E = filtered_train_E, filtered_train_X)
    filtered_var_name = colnames(filtered_train_X)
    rss_vals = numeric(length(filtered_var_name))
    model_data = expanded_train_data_with_weight[valid_indices,]
    for (vi in seq_along(var_name)) {
      variable_name = var_name[vi]
      form = as.formula(paste("E ~", variable_name))
      fit_uni = earth(form,
                      data = model_data,
                      glm = list(family = binomial),
                      weights = filtered_weights)
      rss_vals[vi] = fit_uni$rss
    }
    sig_var_name = filtered_var_name[rank(rss_vals) <= var_threshold]
    gam.model <- as.formula(paste("observed_time~", paste("s(",
                                                          sig_var_name,
                                                          ", k = ", k, ")", sep = "", collapse = "+")))
  }
  if (sum(!cts.x) == length(cts.x)) {
    gam.model <- as.formula(paste("observed_time~", paste(colnames(train_X),
                                                          collapse = "+"), sep = ""))
  }
  cox_gam_model <- mgcv::gam(gam.model, family=mgcv::cox.ph(), data = train_data, weights=sigma, method = "REML")
  test_X_with_time = data.frame(observed_time = time_point, test_X)
  coxph_gam_survival_probability <- predict(cox_gam_model, newdata = test_X_with_time, type = "response")
  coxph_gam_EP = data.frame(Additive_CoxPH = 1-coxph_gam_survival_probability)
  return(coxph_gam_EP)
}



