## ---------------------------------------------------------------------------------------------------------------
true_survival_function <- function(dist, time, params){
  # This function computes the BASELINE true survival probability for different survival distributions at a given time.
  # The supported distributions are exponential, Weibull, log-normal, and log-logistic.
  # Arguments:
  #   dist: A character string specifying the survival distribution to use. 
  #         Options of dist: c('exponential', 'Weibull', 'log-normal', or 'log-logistic')
  #   time: A numeric value representing the time point at which to compute the survival probability.
  #   params: A list of distribution-specific parameters.
  #           - For 'exponential', include 'lambda_base' (rate parameter).
  #           - For 'Weibull', include 'lambda_base' (1/(scale parameter)) and 'alpha' (shape parameter).
  #           - For 'log-normal', include 'mean' (mean of log-time) and 'sd' (standard deviation of log-time).
  #           - For 'log-logistic', include 'lambda_base' (1/(scale parameter)) and 'alpha' (shape parameter).
  #
  # Returns:
  #   A numeric value representing the survival probability at the given time.
  # Below is the formula of baseline survival function
  # Exponential S(t) = exp(-lambda*t)
  # Weibull S(t) = exp(-(lambda*t)^alpha)
  # Log Normal S(t) = 1-CDF((ln(t) - mean)/sd), CDF() here is the CDF of standard normal
  # Log Logistic S(t) = 1/(1+(lambda*t)^alpha)
    if (dist == "exponential") {
    return(exp(-params$lambda_base * time))
  } else if (dist == "Weibull") {
    return(exp(-(params$lambda_base * time)^params$alpha))
  } else if (dist == "log-normal") {
    return(plnorm(time, meanlog = params$mean, sdlog = params$sd, lower = FALSE))
  } else if (dist == "log-logistic") {
    return(1 / (1 +(params$lambda_base * time)^params$alpha))
  } else {
    stop("Error: Unsupported distribution.")
  }
}