## ---------------------------------------------------------------------------------------------------------------
# A helper function, using Kaplan Meier estimator of survival distribution of the censoring times.
Get_Ghat <- function(observed_time, sigma){
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

# Get_Ghat <-function(observed_time, sigma){
#   # A helper function, using Kaplan Meier estimator of survival distribution of the censoring times.
#   km_censor = survfit(Surv(observed_time, 1-sigma)~1, data = full_dt)
#   step_function = stepfun(km_censor$time, c(1, km_censor$surv))
#   return(step_function)
# }

Ghat_newtime <- function(step_function, new_observed_time, time_point){
  ordered_weights = step_function(new_observed_time)
  new_G = ordered_weights[rank(new_observed_time)]
  new_G[new_observed_time >= time_point] = step_function(time_point)
  return(new_G)
}


## ---------------------------------------------------------------------------------------------------------------
#' Compute Inverse Probability of Censoring Weights (IPCW)
#'
#' This function calculates inverse probability of censoring weights (IPCW) 
#' using the Kaplan-Meier estimator to adjust for censoring in survival analysis.
#'
#' @param observed_time Numeric vector of observed event or censoring times.
#' @param sigma Binary censoring indicator (1 = event, 0 = censored).
#' @param time_point The specific cut-off time
#' @param return_type Character string indicating the output format. 
#'   - `"data"` (default): Returns a data frame containing columns observed_time, sigma, and IPCW.
#'   - `"IPCW"`: Returns only the IPCW vector.

get_IPCW_relevant <- function(observed_time, sigma, time_point, return_type = "data"){
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

# get_IPCW_relevant <- function(dt, time_point, return_type = "data"){
#   km_censor_Xi = survfit(Surv(observed_time, 1-sigma)~1, data = dt)
#   survest_Xi = stepfun(km_censor_Xi$time, c(1, km_censor_Xi$surv))
#   censor_prob_Xi = survest_Xi(ifelse(dt$observed_time<time_point, dt$observed_time,time_point))
#   dt$G_hat_Vi = censor_prob_Xi
#   IPCW = ifelse(pmin(dt$event_time, time_point)<dt$censor_time, 1/pmax(dt$G_hat_Vi, 1e-15), 0)
#   if (return_type == "IPCW"){
#     return(IPCW)
#   }
#   else{
#     dt$IPCW = IPCW
#     return(dt)
#   }
# }

## ---------------------------------------------------------------------------------------------------------------
#' Compute Inverse Probability of Censoring Weights (IPCW) with discrete-time modelling
#'
#' This function calculates inverse probability of censoring weights (IPCW) 
#' using the Kaplan-Meier estimator to adjust for censoring in survival analysis and discrete-time modelling
#'
#' @param interval A numeric vector that divides the range from 0 to time_point into equal intervals
#' @param data Dataframe containing covariates, observed_time, sigma and E (status)
#' @param time_point The specific cut-off time
#' @param var_name Names of the covaraites
#' @param return_type Character string indicating the output format. 
#'   - `"model"` (default): Returns a glm model using IPC weighted observations.
#'   - `"dataset"`: Returns a data frame containing columns covariates, observed_time, sigma, E, and IPCW.

bin_combined_ipcw <- function(interval, data, time_point, var_name, return_type="model"){
  total=data.frame()
  data$id = c(1: nrow(data))
  Ghat = Get_Ghat(observed_time = data$observed_time, sigma = data$sigma)
  for(i in 1:(length(interval)-1)){
    start = interval[i]
    end = interval[i+1]
    numerator = Ghat(start)
    batch = data[data$observed_time>start,]
    batch$E = get_status(observed_time = batch$observed_time,
                         sigma = batch$sigma,
                         time_point = end)
    batch$G = Ghat_newtime(step_function = Ghat, new_observed_time = batch$observed_time, time_point = end)
    batch$IPCW = ifelse(is.na(batch$E), 0, numerator/batch$G)
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