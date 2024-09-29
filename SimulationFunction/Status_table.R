## ------------------------------------------------------------------------------------------
# knitr::purl("~/IPCweighting/SimulationFunction/IPCWmethod.Rmd", output = "~/IPCweighting/SimulationFunction/IPCWmethod.R")
# source("~/IPCweighting/SimulationFunction/IPCWmethod.R")


## ------------------------------------------------------------------------------------------
check_time_vs_status <- function(dt,time_point){
  # return a table shows the number of E=0, E=1, E unknown, respectively
  num_E_1 = sum((get_status(dt$observed_time, dt$sigma, time_point)) == 1, na.rm = TRUE)
  num_E_0 = sum((get_status(dt$observed_time, dt$sigma, time_point)) == 0, na.rm = TRUE)
  num_E_unknown = sum(is.na(get_status(dt$observed_time, dt$sigma, time_point)))
  status_table = data.frame(status = c("E = 1", "E = 0", "E unknown"),
                            Count = c(num_E_1, num_E_0, num_E_unknown),
                            stringsAsFactors = FALSE)
  colnames(status_table)[1] = paste("Status_before_year_", time_point, sep = "")
  return(status_table)
}


## ------------------------------------------------------------------------------------------
check_censor_in_time_interval<- function(dt,time_point){
  # return a table shows the number of observations censored between i-1 and i, where i=1,2,..., time_point
  result_table =data.frame(Time_Interval = character(),
                         Count = integer(),
                         stringsAsFactors = FALSE)
  for (i in 1:time_point) {
    count = nrow(dt[dt$sigma == 0 & dt$observed_time >= i - 1 & dt$observed_time <= i, ])
    result_table = rbind(result_table, 
                          data.frame(Time_Interval = paste(i-1, "to", i), 
                                     Count_Censored_Data_in_Time_Tnterval = count))
  }
  return(result_table)
}

