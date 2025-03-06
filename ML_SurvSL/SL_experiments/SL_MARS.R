# library(foreach)
# library(doParallel)
source("~/IPCweighting/ML_SurvSL/sources.R")

ncores = detectCores()
registerDoParallel(cores=ncores-3)# Shows the number of Parallel Workers to be used
print(ncores) # this how many cores are available, and how many you have requested.
getDoParWorkers()# you can compare with the number of actual workers

dist <- "exponential"
params <- list()
params$lambda_base = 0.15
censor_params <- list()
censor_dist <- "uniform"
censor_params$start <- 1
censor_params$end <- 8
num_obs = 250
train_prop = 0.7
neighbors_range = (1:floor(sqrt(num_obs * train_prop)))[(1:floor(sqrt(num_obs * train_prop))) %% 2 == 1]
time_point = 5
num_covariates = 5
measure = "Brier_Score"

set.seed(2003)
X <- matrix(rnorm(num_obs * num_covariates), nrow = num_obs, ncol = num_covariates)
colnames(X) <- paste0("X", 1:num_covariates)
lc = 2 * sin(pi * X[,1]) + 1.3 * (X[,2]^3) - 1.5 * X[,3] * X[,4] - 0.3 * abs(X[,5]) + 0.6 * X[, 4]
true_surv = (true_survival_function(dist, time_point, params = params))^exp(lc)
meta_learner_params = list(penalized = TRUE,
                           alpha_num = 1,
                           lambda_values = 10^seq(0, -5, length = 60),
                           useMin = TRUE)

GlobalFunctions = ls(globalenv())

start_iter = 1
end_iter = 30
max_bin = 10

results_total_list <- foreach(w = start_iter:end_iter, .packages = c("MASS", "dplyr","glmnet", "survival", "caret", "rpart", "earth"), .export=GlobalFunctions) %dopar% {
  set.seed(20+w)
  bin_mars <- lapply(seq(max_bin), function(bin) {
    print(bin)
    learners_glm = create.Learner("dIPCWSL.glm", tune = list(num_intervals = bin),
                                  detailed_names = TRUE)
    learners_ridge = create.Learner("dIPCWSL.glmnet", tune = list(num_intervals = bin, alpha = c(0), lambda_grid = 10^seq(-1, -4, length = 50)),
                                    detailed_names = TRUE)
    learners_elastic = create.Learner("dIPCWSL.glmnet", tune = list(num_intervals = bin, alpha = c(0.5), lambda_grid = 10^seq(-1, -4, length = 50)),
                                      detailed_names = TRUE)
    learners_lasso = create.Learner("dIPCWSL.glmnet", tune = list(num_intervals = bin, alpha = c(1),lambda_grid = 10^seq(-1, -4, length = 50)),
                                    detailed_names = TRUE)
    learners_rpart = create.Learner("dIPCWSL.rpart", tune = list(num_intervals = bin, cp = 2 * 10^seq(-1, -4, length = 40)),
                                    detailed_names = TRUE)
    learners_knn = create.Learner("dIPCWSL.knn", tune = list(num_intervals = bin, k = neighbors_range),
                                  detailed_names = TRUE)
    learners_earth = create.Learner("dIPCWSL.earth", tune = list(num_intervals = bin, degree = c(1), nprune = c(1, 2, 4, 8, 11)),
                                    detailed_names = TRUE)
    SL_library = c(learners_glm$names, learners_ridge$names, learners_elastic$names, learners_lasso$names,
                   learners_rpart$names, learners_knn$names, learners_earth$names)
    run_simulation_SL(num_obs,
                      time_point = time_point,
                      X = X,
                      lc,
                      dist,
                      params,
                      measure = measure,
                      censor_dist,
                      censor_params,
                      train_prop = train_prop,
                      range_intervals = bin,
                      true_surv = true_surv,
                      SL_library = SL_library,
                      meta_learner_params = meta_learner_params,
                      method = "method.glmnet",
                      control = list(),
                      cvControl = list(V = 5, stratifyCV = TRUE, shuffle = TRUE),
                      verbose = FALSE)
  })
  
  results_list = lapply(bin_mars, function(ele) ele$results_df)
  IPCWML_CVRisk = lapply(bin_mars, function(ele) ele$ML_CVRisk)
  IPCWcoef = lapply(bin_mars, function(ele) ele$coef)
  all_layer_testEP = lapply(bin_mars, function(ele) ele$all_layer_testEP)
  if(measure == "C_index"){
    Windex <- which.max(unlist(IPCWML_CVRisk))
  } else {
    Windex <- which.min(unlist(IPCWML_CVRisk))
  }
  c_results = do.call(rbind, results_list)
  Wbest_row = c_results[c_results$Method == paste0("glmnetSLIPCW_", Windex), ]
  Wbest_row$Method = "OPT_glmnetSLIPCW"
  Wbest_row$Opt_bins <- as.numeric(Windex)
  combined_results = rbind(c_results, Wbest_row)
  list(combined_results = combined_results,
       all_layer_testEP=all_layer_testEP,
       IPCWcoef = IPCWcoef)
}

# Combine all results into a single data frame
results_list = lapply(results_total_list, function(item) item$combined_results)
all_layer_testEP = lapply(results_total_list, function(item) item$all_layer_testEP)
IPCWcoef = lapply(results_total_list, function(item) item$IPCWcoef)
combined_df <- bind_rows(results_list)
df <- combined_df %>%
  group_by(Method) %>%
  summarise(across(1:(ncol(combined_df)-1), mean, na.rm = TRUE)) %>%
  arrange(.[[2]]) %>%
  as.data.frame()

saveRDS(list(all_layer_testEP = all_layer_testEP, 
             results_list = results_list, 
             combined_results = df, 
             IPCWcoef = IPCWcoef), paste0("SL_MARS_", dist, "_", num_obs, "_", num_covariates, ".rds"))


