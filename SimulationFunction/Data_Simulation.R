## ----------------------------------------------------------------------------------
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
  else if (dist == "uniform"){
    if (is.null(params$start)&is.null(params$end)){
      stop("Parameters are not enough for uniform function")
    }
  }
  else{
    stop("Unsupported distribution")
  }
}


## ----------------------------------------------------------------------------------
simulation_data <-function(num_obs,
                           event_dist,
                           event_params=list(),
                           censor_dist,
                           censor_params = list(),
                           time_point,
                           X){
  check_event_params = check_params(event_dist, event_params, "event")
  check_censor_params = check_params(censor_dist, censor_params, "censor")
  dt = X
  dt$event_time = sim_T(event_dist, num_obs, params=event_params)
  dt$censor_time = sim_censor_time(censor_dist, num_obs, params = censor_params)
  dt$observed_time = pmin(dt$event_time, dt$censor_time)
  dt$sigma = as.numeric(dt$event_time<dt$censor_time)
  dt$E = get_status(dt$observed_time, dt$sigma, time_point)
  return(dt)
}


## ----------------------------------------------------------------------------------
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


