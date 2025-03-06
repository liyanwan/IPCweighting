# library(foreach)
# library(doParallel)
source("~/IPCweighting/ML_SurvSL/sources.R")

ncores = detectCores()
registerDoParallel(cores=ncores-3)# Shows the number of Parallel Workers to be used
print(ncores) # this how many cores are available, and how many you have requested.
getDoParWorkers()# you can compare with the number of actual workers

# Define parameters
dist <- "log-normal"
params <- list()
censor_params <- list()
censor_dist <- "uniform"
censor_params$start <- 1
censor_params$end <- 8
params$mean <- 1.5
params$sd <- 0.8
time_point = 5

num_obs <- 250
train_prop = 0.7
neighbors_range = (1:floor(sqrt(num_obs * train_prop)))[(1:floor(sqrt(num_obs * train_prop))) %% 2 == 1]
meta_learner_params = list(penalized = TRUE,
                           alpha_num = 1,
                           lambda_values = 10^seq(0, -5, length = 60),
                           useMin = TRUE)
num_covariates <- 100
num_beta = num_covariates
set.seed(2003)
X <- matrix(rnorm(num_obs * num_covariates), nrow = num_obs, ncol = num_covariates)
colnames(X) <- paste0("X", 1:num_covariates)

GlobalFunctions = ls(globalenv())

start_iter = 1
end_iter = 1500
max_bin = 10

results_total_list <- foreach(w = start_iter:end_iter, .packages = c("MASS", "dplyr","glmnet", "survival", "caret", "rpart", "earth"), .export=GlobalFunctions) %dopar% {
  set.seed(20 + w)
  beta <- c(rnorm(num_beta, mean = 0, sd = 0.3), rep(0, num_covariates-num_beta))
  lc = rowSums(sweep(X,2,beta,"*"))
  true_surv = (true_survival_function(dist, time_point, params = params))^exp(lc)
  bin_hdall <- lapply(seq(max_bin), function(bin) {
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
    learners_earth_0.02 = create.Learner("dIPCWSL.earth", tune = list(num_intervals = bin, degree = c(1), nprune = c(2, 4, 6, 8)),
                                         detailed_names = TRUE)
    library_1 = lapply(seq_along(learners_earth_0.02$names), function (m){
      c(learners_earth_0.02$names[m], "screen.corP.s")
    })
    learners_earth_0.1 = create.Learner("dIPCWSL.earth", tune = list(num_intervals = bin, degree = c(1), nprune = c(2,4,6,8,10,12)),
                                        detailed_names = TRUE)
    library_2 = lapply(seq_along(learners_earth_0.1$names), function (m){
      c(learners_earth_0.1$names[m], "screen.corP")
      })
    learners_earth_10 = create.Learner("dIPCWSL.earth", tune = list(num_intervals = bin, degree = c(1), nprune = c(2,4,6,8,10)),
                                       detailed_names = TRUE)
    library_3 = lapply(seq_along(learners_earth_10$names), function (m){
      c(learners_earth_10$names[m], "screen.corRank")
      })
    learners_earth_glmnet = create.Learner("dIPCWSL.earth", tune = list(num_intervals = bin, degree = c(1), nprune = c(4, 8, 10, 12, 14, 16, 20)),
                                           detailed_names = TRUE)
    library_4 = lapply(seq_along(learners_earth_glmnet$names), function (m){
      c(learners_earth_glmnet$names[m], "screen.glmnet")
      })
    SL_library = c(as.list(learners_glm$names),
                   as.list(learners_ridge$names),
                   as.list(learners_elastic$names),
                   as.list(learners_rpart$names),
                   as.list(learners_knn$names),
                   library_1,
                   library_2,
                   library_3,
                   library_4)
    
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
  
  results_list = lapply(bin_hdall, function(ele) ele$results_df)
  IPCWML_CVRisk = lapply(bin_hdall, function(ele) ele$ML_CVRisk)
  IPCWcoef = lapply(bin_hdall, function(ele) ele$coef)
  all_layer_testEP = lapply(bin_hdall, function(ele) ele$all_layer_testEP)
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
             IPCWcoef = IPCWcoef), paste0("SL_hdall_", dist, "_", num_obs, "_", num_covariates, ".rds"))


