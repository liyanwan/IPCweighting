## ---------------------------------------------------------------------------------------------------------------
ols_error <- function(true_surv, est_surv) {
  return(mean((true_surv - est_surv)^2))
}


## ---------------------------------------------------------------------------------------------------------------
weighted_loglikelihood <- function(IPCW, status, event_prob) {
  event_prob = pmax(pmin(event_prob, 1 - 1e-6), 1e-6)
  log_likelihood = -mean(
    ifelse(IPCW == 0,
           0,
           IPCW*(status * log(event_prob) + (1-status)*log(1 - event_prob))
    )
  )
  return(log_likelihood)
}


## ---------------------------------------------------------------------------------------------------------------
Weighted_Brier_Score <- function(event_prob, status, IPCW){
  score = mean(ifelse(IPCW==0,
                      0,
                      IPCW * (event_prob - status)^2))
  return(score)
}
