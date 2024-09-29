## -----------------------------------------------------------------------------------------------------
rm(list = ls())
library(survival)
library(dplyr)
library(caret)
library(ggplot2)
library(knitr)
knitr::purl("~/IPCweighting/SimulationFunction/IPCWmethod.Rmd", output = "~/IPCweighting/SimulationFunction/IPCWmethod.R")
source("~/IPCweighting/SimulationFunction/IPCWmethod.R")
knitr::purl("~/IPCweighting/IPCW/estimation_naive_cox_ipcw.Rmd", output = "~/IPCweighting/IPCW/estimation_naive_cox_ipcw.R")
source("~/IPCweighting/IPCW/estimation_naive_cox_ipcw.R")


## -----------------------------------------------------------------------------------------------------
bin_combined_ipcw<-function(interval,train_data){
  # This is used to get the model of "binned" IPCW
  total=data.frame()
  for(i in 1:(length(interval)-1)){
    start = interval[i]
    end = interval[i+1]
    numerator = Ghat(start)
    batch = train_data[train_data$observed_time>start,]
    batch$E = NULL
    batch$E = get_status(batch$observed_time,batch$sigma,end)
    batch$G = ifelse(batch$observed_time<end, Ghat(batch$observed_time),Ghat(end))
    batch$IPCW = ifelse(pmin(batch$event_time, end)<batch$censor_time, 
                                  numerator/batch$G, 
                                  0)
    if (ncol(total) == 0) {
      total = data.frame(batch)
    } else {
      total = rbind(total, batch)
    }
  }
  formula_str = paste("E ~", paste(var_name, collapse = " + "))
  formula = as.formula(formula_str)
  model = glm(formula, data = total, family = binomial, weights = IPCW)
  return(model)
}


## -----------------------------------------------------------------------------------------------------
bin_surv_probability <-function(model_list=list(), data, var_name){
  pred_survival_probability = data.frame()
  for(i in 1:length(model_list)){
    name_of_bin = paste0("Bin_", i)
    # survprob_in_bin = 1 - predict(model, newdata = data[var_name], type = "response")
    survprob_in_bin = data.frame(setNames((1 - predict(model_list[[i]], newdata = data[var_name], type = "response")), name_of_bin))
    if (ncol(pred_survival_probability) == 0) {
      pred_survival_probability = survprob_in_bin
    }
    else{
      pred_survival_probability = cbind(pred_survival_probability, survprob_in_bin)
    }
  }
  return(pred_survival_probability)
}


## -----------------------------------------------------------------------------------------------------
bin_separate_ipcw <-function(interval,train_data){
  # This is used to get the model of "binned" IPCW
  model_list = list()
  for(i in 1:(length(interval)-1)){
    start = interval[i]
    end = interval[i+1]
    numerator = Ghat(start)
    batch = train_data[train_data$observed_time>start,]
    batch$E = NULL
    batch$E = get_status(batch$observed_time,batch$sigma,end)
    batch$G = ifelse(batch$observed_time<end, Ghat(batch$observed_time),Ghat(end))
    batch$IPCW = ifelse(pmin(batch$event_time, end)<batch$censor_time, 
                        numerator/batch$G, 
                        0)
    formula = as.formula(paste("E ~", paste(var_name, collapse = " + ")))
    batch_model = glm(formula, data = batch, family = binomial, weights = IPCW)
    model_list[[i]] = batch_model
  }
  return(model_list)
}

