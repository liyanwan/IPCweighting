## -----------------------------------------------------------------------------------------------------
rm(list = ls())
library(survival)
library(dplyr)
library(caret)
library(ggplot2)
library(knitr)
knitr::purl("~/IPCweighting/IPCWmethod.Rmd", output = "~/IPCweighting/IPCWmethod.R")
source("~/IPCweighting/IPCWmethod.R")
knitr::purl("~/IPCweighting/estimation_naive_cox_ipcw.Rmd", output = "~/IPCweighting/estimation_naive_cox_ipcw.R")
source("~/IPCweighting/estimation_naive_cox_ipcw.R")


## -----------------------------------------------------------------------------------------------------
num_obs = 1000
dist = "exponential"
params = list()
time_point = 5
if(dist=="exponential"){
  params$lambda_base = 0.5
  lambda_base = 0
}else if(dist == "Weibull"){
  params$lambda_base = 2
  params$alpha = 0.5
} else if(dist == "log-logistic"){
  params$lambda_base = 4
  params$alpha = 1.2
} else if (dist == "log-normal"){
  params$mean = 2
  params$sd = 0.7
} else{
  stop("Unsupported distribution")
}


## -----------------------------------------------------------------------------------------------------
# Simulate covariates and event_time
X = matrix(ncol = 0, nrow = 0)
mean_X = c()
X = get_cont_X(X, 'uniform',num_obs,1,4,"X1")
X = get_cont_X(X, 'uniform',num_obs,0.2,1,"X2")
beta = sim_beta(0.6, time_point, mean_X) #sim_beta("uniform", length(mean_X),list(start = -2, end = 2))
if (ncol(X)!=length(beta)){
  stop("X and beta are not multiplicable")
}
params$lc = rowSums(sweep(X,2,beta,"*"))
power = exp(rowSums(sweep(X,2,beta,"*")))
dt = data.frame(power = power)
dt = cbind(dt,X)
var_name = colnames(X)
dt$event_time = sim_T(dist, num_obs, params=params, TRUE)


## -----------------------------------------------------------------------------------------------------
# simulate censoring time
censor_params = list()
censor_dist = "uniform"
if(censor_dist=="exponential"){
  censor_params$lambda = 0.5
}else if(censor_dist == "Weibull"){
  censor_params$lambda = 2
  censor_params$alpha = 0.5
} else if(censor_dist == "log-logistic"){
  censor_params$lambda = 4
  censor_params$alpha = 1.2
} else if (censor_dist == "log-normal"){
  censor_params$mean = 2
  censor_params$sd = 0.7
} else if (censor_dist == "uniform"){
  censor_params$start = 1
  censor_params$end = 7
} else{
  stop("Unsupported distribution")
}
dt$censor_time = sim_censor_time(censor_dist, num_obs, params = list(start = 2.5, end = 9))


## -----------------------------------------------------------------------------------------------------
dt$observed_time = pmin(dt$event_time, dt$censor_time)
dt$sigma = as.numeric(dt$event_time<dt$censor_time)
dt$E = get_status(dt$observed_time, dt$sigma, time_point)
print(check_time_vs_status(dt,time_point))
print(check_censor_in_time_interval(dt,time_point))
result_table =data.frame(Time_Interval = character(),
                       Count = integer(),
                       stringsAsFactors = FALSE)
for (i in 1:time_point) {
  count = nrow(dt[dt$observed_time >= i - 1 & dt$observed_time <= i, ])
  result_table = rbind(result_table, 
                        data.frame(Time_Interval = paste(i-1, "to", i), 
                                   Number_of_Observations_in_Time_Tnterval = count))
}
print(result_table)


## -----------------------------------------------------------------------------------------------------
multilayer_ipcw<-function(interval,train_data){
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
Ghat = Get_Ghat(time_point, dt)
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


## ----warning=FALSE------------------------------------------------------------------------------------
# Cross Validation, binning intervals equally
max_intervals = 30
k = 5
Ghat = Get_Ghat(time_point, dt)
folds = createFolds(dt$event_time, k = k)
results = data.frame(intervals = 1:max_intervals, performance_combined = NA, performance_separate = NA)
for (num_intervals in 1:max_intervals) {
  fold_performance_separate = numeric(k)
  fold_performance_combined = numeric(k)
  for (f in 1:k) {
    train_data = dt[-folds[[f]], ]
    test_data = dt[folds[[f]], ]
    interval = seq(0, time_point, length.out = num_intervals + 1)
    model_list = bin_separate_ipcw(interval, train_data)
    pred_survival_probability = bin_surv_probability(model_list, test_data, var_name)
    test_data$pred_separate_prob = apply(pred_survival_probability, 1, prod)
    test_true_surv = (true_survival_function(dist, time_point, params = params))^exp(rowSums(sweep(test_data[var_name],2,beta,"*")))
    model = multilayer_ipcw(interval, train_data)
    test_data$pred_combine_prob = 1 - predict(model, newdata = test_data[var_name], type = "response")
    fold_performance_combined[f] = ols_error(test_true_surv, (test_data$pred_combine_prob)^num_intervals)
    fold_performance_separate[f] = ols_error(test_true_surv, (test_data$pred_separate_prob))
}
  results$performance_combined[num_intervals] = mean(fold_performance_combined)
  results$performance_separate[num_intervals] = mean(fold_performance_separate)
  
}
best_intervals_combined = results$intervals[which.min(results$performance_combined)]
best_intervals_separate = results$intervals[which.min(results$performance_separate)]

print(paste("Optimal number of intervals combined:", best_intervals_combined))
print(paste("Optimal number of intervals separate:", best_intervals_separate))


## -----------------------------------------------------------------------------------------------------
# Plot of performance with respect to the number of bins
plot(x=c(1:max_intervals), y = results$performance_combined, xlab = "Number of Intervals", ylab = "OLS_Error")
plot(x=c(1:max_intervals), y = results$performance_separate, xlab = "Number of Intervals", ylab = "OLS_Error")


## -----------------------------------------------------------------------------------------------------
# Compute 1 bin IPCW ols_error
true_surv = (true_survival_function(dist, time_point, params = params))^power
ipcw_est_surv = ipcw_estimate(dt,dt, time_point, var_name)
ols_ipcw = ols_error(true_surv,ipcw_est_surv)
ols_ipcw

# Compute optimal_bin IPCW ols_error
opt_interval = seq(0, time_point, length.out = best_intervals_combined + 1)
multi_ipcw_model = multilayer_ipcw(opt_interval, dt)
multi_layer_ipcw_pred_prob = 1-(predict(multi_ipcw_model, newdata = dt[var_name], type = "response"))
ols_bin_combine_ipcw = ols_error(true_surv, (multi_layer_ipcw_pred_prob)^best_intervals_combined)
ols_bin_combine_ipcw

opt_interval_separate = seq(0, time_point, length.out = best_intervals_separate + 1)
model_list = bin_separate_ipcw(opt_interval_separate, dt)
pred_survival_probability = bin_surv_probability(model_list, dt, var_name)
pred_separate_prob = apply(pred_survival_probability, 1, prod)
ols_bin_separate_ipcw = ols_error(true_surv, pred_separate_prob)
ols_bin_separate_ipcw

