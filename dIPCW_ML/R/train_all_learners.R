## ---------------------------------------------------------------------------------------------------------------
#' Build a weighted model and estimate event probabilities on X_validation using cross-validation
#' results for each learner function.
#'
#' @param train_X A dataframe of covariates.
#' @param train_Y A dataframe containing observed_time, sigma, and E.
#' @param X_validation A dataframe of covariates used for validation.
#' @param time_point The specified cut-off time for evaluation.
#' @param measure The metric used for selecting the optimal tuning parameters.
#' @param range_intervals A numeric vector defining the number of intervals to evaluate during cross-validation.
#' @param learner_list A list of learner functions to be used in the model.
#' @param params_list A list of parameter sets corresponding to each learner function.
#' @param foldid A vector specifying the fold assignment for each observation.
#' @param ... Additional parameters that can be passed to each learner function.

train_all_learners_all <- function(train_X, train_Y, X_validation = NULL, time_point, measure, range_intervals, learner_list,params_list, foldid = NULL, ...) {
  model_list_with_bin <- vector("list", length(learner_list))
  for (j in seq_along(learner_list)){
    learner = learner_list[[j]]
    parms = params_list[[j]]
    learner_info = do.call(learner, c(list(X = train_X, Y = train_Y, 
                                           X_validation = X_validation, 
                                           range_intervals = range_intervals, 
                                           measure = measure, foldid = foldid,
                                           time_point = time_point), ..., parms)) # ALL params except proxy are given
    model_list_with_bin[[j]] <- list(model = learner_info$model,
                                     bin = learner_info$optbin,
                                     test_event_pred = learner_info$pred,
                                     cv = learner_info$results)
    model_list_with_bin[[j]]$class_name = learner_info$class_name
  }
  return(model_list_with_bin)
}


## ---------------------------------------------------------------------------------------------------------------
#' Build a weighted model and estimate event probabilities on X_validation using cross-validation
#' results for each learner function with 1-binned interval.

fit_singlebin <- function(train_X, train_Y, X_validation, time_point, measure, learner_list, params_list, return_type = "single", foldid = NULL, ...){
  train_data = cbind(train_X, train_Y)
  single_test_event_probability = list()
  model_list_with_single_bin = train_all_learners_all(
    train_X = train_X, train_Y = train_Y, X_validation = X_validation, learner_list, params_list,
    time_point = time_point, range_intervals = c(1), foldid = foldid, measure = measure, ...
  )
  single_model_event_probability_layer = data.frame()
  for (i in seq_along(model_list_with_single_bin)) {
    pred_test_event_probability = as.vector(model_list_with_single_bin[[i]]$test_event_pred)
    class_name = model_list_with_single_bin[[i]]$class_name
    col_name = if (return_type == "naive") gsub("Single", "Naive", class_name) else class_name
    if (ncol(single_model_event_probability_layer) == 0) {
      single_model_event_probability_layer = data.frame(pred_test_event_probability)
      colnames(single_model_event_probability_layer) = col_name
    } else if (col_name %in% colnames(single_model_event_probability_layer)) {
      newname = paste0(col_name,'.', i)
      single_model_event_probability_layer[[newname]] = pred_test_event_probability
    } else {
      single_model_event_probability_layer[[col_name]] = pred_test_event_probability
    }
  }
  return(list(single_test_event_probability = single_model_event_probability_layer))
}