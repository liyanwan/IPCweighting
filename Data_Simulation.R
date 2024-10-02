## --------------------------------------------------------------------------------------------------
knitr::purl("~/IPCweighting/IPCWmethod.Rmd", output = "~/IPCweighting/IPCWmethod.R")
source("~/IPCweighting/IPCWmethod.R")


## --------------------------------------------------------------------------------------------------
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


## --------------------------------------------------------------------------------------------------
Simulation_data <-function(num_obs,
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

