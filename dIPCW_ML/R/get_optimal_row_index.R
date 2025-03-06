#' This function identifies optimal row indices from a dataset based on specific criteria.
#' The function works with datasets containing 4, 5, or more columns and selects rows
#' based on the highest or lowest values of certain columns. It also provides an option
#' to use a standard error (SE) threshold instead of exact min/max values.

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





