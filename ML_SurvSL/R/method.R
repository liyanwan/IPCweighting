# outline for SuperLearner methods
# these should always have class 'SL.method'
#
# The SuperLearner method is a coupling of the estimation algorithm for the algorithm weights (coefficients) and the model to combine the algorithms
#
# 2 parts need to be included:
#   1) compute coefficients
#   2) compute predictions

method.template <- function() {
  out <- list(
  # require allows you to pass a character vector with required packages
  # use NULL if no required packages
  require = NULL,
  # computeCoef is a function that returns a list with three elements:
  # 1) coef: the weights (coefficients) for each algorithm
  # 2) cvRisk: the V-fold CV risk for each algorithm
  # 3) optimizer: (optional) the result object from the optimization of the weights.
  computeCoef = function(Z, Y, libraryNames, trainobsWeights, control, verbose, ...) {
    cvRisk <- numeric()
    coef <- numeric()
    out <- list(cvRisk = cvRisk, coef = coef, optimizer = NULL)
    return(out)
  },
  # computePred is a function that takes the weights and the predicted values from each algorithm in the library and combines them based on the model to output the super learner predicted values
  computePred = function(predE, coef, control, ...) {
    out <- crossprod(t(predE), coef)
    return(out)
  }
  )
  invisible(out)
}

write.method.template <- function(file = '', ...) {
  cat('method.template <- function() {\n  out <- list(\n    # require allows you to pass a character vector with required packages\n    # use NULL if no required packages\n    require = NULL,\n\n    # computeCoef is a function that returns a list with two elements:\n    # 1) coef: the weights (coefficients) for each algorithm\n    # 2) cvRisk: the V-fold CV risk for each algorithm\n    computeCoef = function(Z, E, libraryNames, trainobsWeights, control, verbose, ...) {\n      cvRisk <- numeric()\n      coef <- numeric()\n      out <- list(cvRisk = cvRisk, coef = coef)\n      return(out)\n    },\n\n    # computePred is a function that takes the weights and the predicted values\n    # from each algorithm in the library and combines them based on the model to\n    # output the super learner predicted values\n    computePred = function(predE, coef, control, ...) {\n      out <- crossprod(t(predE), coef)\n      return(out)\n    }\n    )\n    invisible(out)\n  }', file = file, ...)
}


method.glmnet <- function() {
  computePred = function(predE, coef, ...) {
    if(length(coef) > 1){
      intercept = coef[1]
      coef = coef[-1]
      if (sum(coef != 0) == 0) {
        warning("All metalearner coefficients are zero, predictions will all be 0", call. = FALSE)
      }
      plogis(intercept + (as.matrix(predE[, coef != 0]) %*%
                          matrix(coef[coef != 0])))
    }
  }
  computeCoef = function(Z, Y, libraryNames, trainobsWeights, alpha_num, 
                         lambda_grid, folds, time_point, useMin, measure = "C-index",...) {
    # check for duplicated columns
    # set a tolerance
    E = Y$E
    noNA_indices = which(!is.na(E))
    E_noNA = E[noNA_indices]
    Z_noNA = Z[noNA_indices, , drop = FALSE]
    if(dim(Z_noNA)[2] == 1){
      stop("method.glmnet requires at least two candidates in SL.library")
    }
    trainobsWeights_noNA = trainobsWeights[noNA_indices]
    intercept_list = list()
    coef_list = list()
    meta_learner_list = list()
    tol = 8
    dupCols <- which(duplicated(round(Z_noNA, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    if(anyDupCols){
      # if present, throw warning identifying learners
      warning(paste0(paste0(libraryNames[dupCols],collapse = ", "),
                     " are duplicates of previous learners.",
                     " Removing from super learner."))
    }
    cvRisk <- apply(Z_noNA, 2, function(x) -sum(2 * trainobsWeights_noNA *
                                                  ifelse(E_noNA, plogis(x, log.p=TRUE),
                                                         plogis(x, log.p=TRUE, lower.tail=FALSE))))
    names(cvRisk) <- libraryNames
    if (is.null(alpha_num)){
      meta_learner = glm(y = E, x = Z, family = binomial, weights = trainobsWeights)
      intercept = as.vector(coefficients(meta_learner))[1]
      coef = as.vector(coefficients(meta_learner))[-1]
      if (anyNA(coef)) {
        warning("Some algorithms have weights of NA, setting to 0.")
        coef[is.na(coef)] <- 0
      }
      intercept_list = list(intercept_list)
      coef_list = list(coef)
      meta_learner_list = list(meta_learner)
    } else{
      if (is.matrix(Z)) {
        Z = as.data.frame(Z)
      }
      dt = data.frame(Z, Y)
      var_name = colnames(Z)
      results =  cross_validate_lambda_SingleMeasure(X = Z,
                                                     Y = Y,
                                                     time_point,
                                                     alpha_num = alpha_num,
                                                     lambda_values = lambda_grid,
                                                     range_intervals = c(1), 
                                                     k = max(folds), 
                                                     foldid = folds,
                                                     proxy_data = NULL,
                                                     measure = measure)
      optimal_results = results$optimal_results
      optimal_index_measure_rows = get_optimal_row_indices_SingleMeasure(optimal_results, useMin = useMin, measure = measure)
      optimal_row = optimal_index_measure_rows[[1]]
      ML_CVRisk = optimal_results[optimal_row, 3]
      if(length(optimal_row) != 1){
        stop("GLMNET error: optimal index should be length 1")
      }
      optimal_lambda = optimal_results[optimal_row, 2]
      meta_learner <- glmnet(x = Z_noNA, 
                             y = E_noNA, 
                             family = "binomial", 
                             weights = trainobsWeights_noNA, 
                             alpha = alpha_num, 
                             lambda = lambda_grid)
      intercept = as.vector(coef(meta_learner, s = optimal_lambda)[1])
      coef = as.vector(coef(meta_learner, s = optimal_lambda)[-1])
      if (anyNA(coef)) {
        warning("Some algorithms have weights of NA, setting to 0.")
        coef[is.na(coef)] <- 0
      }
      }
    out <- list(cvRisk = cvRisk, intercept = intercept, coef = coef, 
                meta_learner = meta_learner, optimal_results = optimal_results,
                ML_CVRisk = ML_CVRisk)
    return(out)
    }
  
  list(require = "glmnet",
       computeCoef = computeCoef,
       computePred = computePred)
}
