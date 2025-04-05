source("~/IPCweighting/dIPCW_ML/sources.R")
source("~/IPCweighting/dIPCW_ML/library.R")

# ncores = 25
# registerDoParallel(cores=ncores)# Shows the number of Parallel Workers to be used
# print(ncores) # this how many cores are available, and how many you have requested.
# getDoParWorkers()# you can compare with the number of actual workers

# Define parameters
dist <- "log-normal"
params <- list()
params$mean <- 1.5
params$sd <- 0.8
#params$lambda_base = 0.15
censor_params <- list()
censor_dist <- "uniform"
censor_params$start <- 1
censor_params$end <- 8
learner_list = list(base_mars_BinnedIPCW_all)
params_list <- list(
  list(nprune_values = seq(2, 20, by = 2), degree_values = c(1, 2), useMin = TRUE)
)

# covariates information
num_obs <- 250
num_covariates <- 25
num_beta <- 25
params$mean <- 1.5
params$sd <- 0.8
#params$lambda_base = 0.15
time_point = 5
train_prop = 0.7
set.seed(2003)
X <- matrix(rnorm(num_obs * num_covariates), nrow = num_obs, ncol = num_covariates)
colnames(X) <- paste0("X", 1:num_covariates)

#defining risk score
lc = 0.2 * sin(pi * X[,1]) + 0.05 * (X[,2]^3) - (0.28 * (X[,3]^2-0.2)) + 0.3 * abs(X[,5]) + 0.15 * X[, 4] # risk score
summary(lc)
true_surv = (true_survival_function(dist, time_point, params = params))^exp(lc)
hist(true_surv)

results_total_list <- lapply(1:1000, function(i) {
  set.seed(20+i)
  params$lc <- lc
  dt <- simulation_data(num_obs, dist, params, censor_dist, censor_params, time_point, X)
  dt$event_time
  if(is.infinite(max(dt$event_time))){
    print(paste0("infinite", i))
    break
  }
})

GlobalFunctions = ls(globalenv())
start_iter = 1
end_iter = 50

results_total_list <- foreach(i = start_iter:end_iter,
                              .packages = c("MASS", "dplyr","glmnet", "survival", "caret", "rpart", "mgcv", "yaImpute", "earth"),
                              .export=GlobalFunctions) %dopar% {
                                dIPCW_ML_split(time_point = time_point,
                                               X = X,
                                               Y = Y,
                                               train_prop = train_prop,
                                               measure = "Brier_Score",
                                               range_intervals = 4,
                                               true_surv = true_surv,
                                               learner_list = learner_list,
                                               params_list = params_list,
                                               proxy_data = NULL,
                                               learner_k = 5,
                                               learner_foldid = NULL,
                                               include_naive = TRUE,
                                               include_SA = TRUE,
                                               include_opt = TRUE,
                                               max_intervals = 10,
                                               min_intervals = 1)
                                }

results_list = lapply(results_total_list, function(item) item$results_df)
all_layer_testEP <- lapply(results_total_list, function(item) item$all_layer_testEP)
combined_df <- bind_rows(results_list)
brier_name = colnames(results_list[[1]])[6]
df <- combined_df %>%
  group_by(Method) %>%
  summarise(across(1:7, \(x) mean(x, na.rm = TRUE))) %>%  # Compute mean for columns 2 to 7
  arrange(.data[[brier_name]]) %>%  # Arrange by the second column
  as.data.frame()



saveRDS(list(results_list = results_list,
             all_layer_testEP = all_layer_testEP,
             combined_results = df), paste0("CM_mars_",dist, "_", num_obs, "_", num_covariates,".rds"))

