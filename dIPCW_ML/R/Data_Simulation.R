## ---------------------------------------------------------------------------------------------------------------
sim_censor_time <- function(dist = "exponential", n, params = list()){
  # Simulate censoring time
  # Arguments:
  #   dist: A character string specifying the distribution to use for simulating the censoring times.
  #         Options include 'exponential', 'Weibull', 'log-logistic', 'log-normal', and 'uniform'.
  #   n: An integer specifying the number of censoring times to generate.
  #   params: A list of distribution-specific parameters:
  #   - For 'exponential', include 'lambda' (rate parameter).
  #   - For 'Weibull', include 'lambda' (1/scale parameter) and 'alpha' (shape parameter).
  #   - For 'log-logistic', include 'lambda' (1/scale parameter) and 'alpha' (shape parameter).
  #   - For 'log-normal', include 'mean' (mean of the log-scale parameter) and 'sd' (standard deviation of the log-scale parameter).
  #   - For 'uniform', include 'start' (minimum value) and 'end' (maximum value).
  # Returns:
  #   A numeric vector of length n representing the simulated censoring times based on the specified distribution.
  if (dist == "exponential") {
    if (!is.null(params$lambda_base)) {
      return(rexp(n, rate = params$lambda_base))
    }
    else{
      stop("Please input the correct parameters.")
    }
  } 
  else if (dist == "Weibull") {
    if (!is.null(params$lambda_base) & !is.null(params$alpha)) {
      return(rweibull(n, shape = params$alpha, scale = 1/params$lambda_base))
    }
    else{
      stop("Please input the correct parameters.")
    }
  } 
  else if (dist == "log-logistic") {
    # Log-logistic distribution using qlogis for quantiles
    if (!is.null(params$lambda_base) & !is.null(params$alpha)) {
      u <- runif(n)
      return((1/(params$lambda_base))* (u / (1 - u))^(1 / params$alpha))
    }
    else{
      stop("Please input the correct parameters.")
    }
  } 
  else if (dist == "log-normal") {
    if (!is.null(params$mean) & !is.null(params$sd)) {
      return(rlnorm(n, meanlog = params$mean, sdlog = params$sd))
    }
    else{
      stop("Please input the correct parameters.")
    }
  } 
  else if (dist == "uniform"){
    if(!is.null(params$start) & !is.null(params$end)){
      censor_time = runif(n, params$start,params$end)
    }
    else{
      stop("Please input the correct parameters.")
    }
  }
  else if (dist == "gamma"){
    if(!is.null(params$shape) & !is.null(params$scale)){
      censor_time = rgamma(n, shape = params$shape, scale = params$scale)
    }
  }
  else {
    stop("Unsupported distribution")
  }
}

## ---------------------------------------------------------------------------------------------------------------
sim_T <- function(dist, n, params=list(), U_notuni){
  # Simulate event time based on baseline survival function.
  # Arguments:
  #   dist: A character string specifying the distribution for simulating event times.
  #         Options are c('exponential', 'Weibull', 'log-logistic', and 'log-normal')
  #   n: An integer specifying the number of event times to simulate.
  #   params: 
  #     - For 'exponential', include 'lambda_base' (rate parameter) and 'lc' (linear combination of covariates).
  #     - For 'Weibull', include 'lambda_base' (1/scale parameter), 'alpha' (shape parameter), and 'lc' (linear combination of covariates).
  #     - For 'log-logistic', include 'lambda_base' (1/scale parameter), 'alpha' (shape parameter), and 'lc' (linear combination of covariates).
  #     - For 'log-normal', include 'mean' (mean of the log-scale), 'sd' (standard deviation of the log-scale), and 'lc' (linear combination of covariates).
  # Returns:
  #   A numeric vector representing the simulated event times based on the specified distribution.
  
  # If baseline survival follows exponential distribution, event time is simulated directly
  # In other cases, simulate events times by suing the inverse of the cumulative hazard function (Using helper function simulate_U)
  if (dist == "exponential") {
    if (!is.null(params$lambda_base) & !is.null(params$lc)) {
      return(rexp(n, rate = (params$lambda_base*exp(params$lc))))
    }
  } else if (dist == "Weibull") {
    if (!is.null(params$lambda_base) & !is.null(params$alpha)& !is.null(params$lc)) {
      U = runif(n, 0, 1)
    }
    return(1/params$lambda_base*(-log(U)/exp(params$lc))^(1/params$alpha))
  }
  else if (dist == "log-logistic") {
    if (!is.null(params$lambda_base) & !is.null(params$alpha)) {
      U = runif(n, 0, 1)
      return((1/(params$lambda_base)) * ((1-U^(1/exp(params$lc))) / U^(1/exp(params$lc)))^(1 / params$alpha))
    }
  } else if (dist == "log-normal") {
    if (!is.null(params$mean) & !is.null(params$sd)& !is.null(params$lc)) {
      U = runif(n, 0, 1)
      return(exp(qnorm(1 - U^(1 / exp(params$lc))) * params$sd + params$mean))
    }
  } else if (dist == "gompertz") {
    if (!is.null(params$shape) & !is.null(params$lambda_base)& !is.null(params$lc)) {
      U = runif(n, 0, 1)
      return(1/params$shape*log(1-(params$shape*log(U))/(params$lambda_base*exp(params$lc))))
    }
  }else {
    stop("Unsupported distribution")
  }
}


## ---------------------------------------------------------------------------------------------------------------
get_status <- function(observed_time, sigma, time_point){
  # Calculate the event status at `time_point`, commly denoted as E
  # Arguments:
  #   observed_time: Numeric vector of observed times.
  #   sigma: Numeric vector indicating censoring status (1: uncensored, 0: censored).
  #   time_point: Numeric value specifying the time point of interest.
  # Returns:
  #   A numeric vector of the same length as `observed_time`, where each element indicates the event status:
  #   - 1 if the event happened prior to time_point.
  #   - 0 if the event did not occur prior to time_point (the observation survived at time_point).
  #   - unknown if the observation is censored before time_point.
  E_obs <- ifelse(observed_time <= time_point & sigma == 1, 1, ifelse(observed_time > time_point, 0, NA))
  return(E_obs)
}


## ------------------------------------------------------------------------------------------------------------
check_params <-function(dist, params=list(), type=c("event","censor")){
  if(dist=="exponential"){
    if (((type == "event") & (is.null(params$lambda_base))) | 
        ((type == "censor") & (is.null(params$lambda)))){
      stop("Parameter lambda_base is needed for exponential function")
    }
  }
  else if(dist == "Weibull"){
    if (((type == "event") & (is.null(params$lambda_base)) & is.null(params$alpha)) | 
        ((type == "censor") & (is.null(params$lambda)) & is.null(params$alpha))){
      stop("Parameters are not enough for Weibull function")
    }
  } 
  else if(dist == "log-logistic"){
    if (((type == "event") & (is.null(params$lambda_base)) & is.null(params$alpha)) | 
        ((type == "censor") & (is.null(params$lambda)) & is.null(params$alpha))){
      stop("Parameters are not enough for log-logistic function")
    }
  } 
  else if (dist == "log-normal"){
    if (is.null(params$mean)&is.null(params$sd)){
      stop("Parameters are not enough for log-normal function")
    }
  }
  else if (dist == "gompertz"){
    if (is.null(params$shape)&is.null(params$lambda_base)){
      stop("Parameters are not enough for log-normal function")
    }
  }
  else if (dist == "uniform"){
    if (is.null(params$start)&is.null(params$end)){
      stop("Parameters are not enough for uniform function")
    }
  }
  else if (dist == "gamma"){
    if (is.null(params$shape)&is.null(params$scale)){
      stop("Parameters are not enough for uniform function")
    }
  }
  else{
    stop("Unsupported distribution")
  }
}


## ------------------------------------------------------------------------------------------------------------
simulation_data <-function(num_obs,
                           event_dist,
                           event_params=list(),
                           censor_dist,
                           censor_params = list(),
                           time_point,
                           X){
  check_event_params = check_params(event_dist, event_params, "event")
  check_censor_params = check_params(censor_dist, censor_params, "censor")
  dt = data.frame(X)
  dt$event_time = sim_T(event_dist, num_obs, params=event_params)
  dt$censor_time = sim_censor_time(censor_dist, num_obs, params = censor_params)
  dt$observed_time = pmin(dt$event_time, dt$censor_time)
  dt$sigma = as.numeric(dt$event_time<dt$censor_time)
  dt$E = get_status(dt$observed_time, dt$sigma, time_point)
  return(dt)
}


## ------------------------------------------------------------------------------------------------------------
plot_pdf <-function(dist, params = list(), time_point){
  # Plot the Probability Density Function (PDF) for Various Distributions
  # Arguments:
  #   dist: A character string specifying the distribution to plot. 
  #         Options are c('exponential', 'Weibull', 'log-logistic', and 'log-normal')
  #   params: A list of distribution-specific parameters:
  #           - For 'exponential', include 'lambda_base' (rate parameter).
  #           - For 'Weibull', include 'lambda_base' (1/scale parameter) and 'alpha' (shape parameter).
  #           - For 'log-logistic', include 'lambda_base' (1/scale parameter) and 'alpha' (shape parameter).
  #           - For 'log-normal', include 'mean' (mean of the log-scale) and 'sd' (standard deviation of the log-scale).
  #   time_point: A numeric value as the upper limit of the x-axis.
  # Returns:
  #   A ggplot object displaying the PDF of the specified distribution.
  x_values = seq(0, time_point, by = 0.01)
  if(dist == "exponential"){
    pdf = dexp(x_values, rate = params$lambda_base)
  } else if(dist == "Weibull"){
    pdf = dweibull(x_values, shape = params$alpha, scale = 1/params$lambda_base)
  } else if (dist == "log-logistic"){
    pdf =(params$lambda_base^params$alpha*params$alpha*x_values^(params$alpha-1))/(1+(params$lambda_base*x_values)^params$alpha)^2
  } else if (dist == "log-normal"){
    pdf = dlnorm(x_values, meanlog = params$mean, sdlog = params$sd)
  }
  plotting = data.frame(x_values = x_values, pdf_values = pdf)
  ggplot(plotting, aes(x = x_values, y = pdf_values)) +
    geom_line(color = "blue", linewidth = 1) +
    labs(title = paste("PDF of", dist), x = "x", y = "Density") +
    theme_minimal()
}


