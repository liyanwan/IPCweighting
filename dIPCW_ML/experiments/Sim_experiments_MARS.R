source("~/IPCweighting/dIPCW_ML/sources.R")
source("~/IPCweighting/dIPCW_ML/library.R")

# ncores = 25
# registerDoParallel(cores=ncores)# Shows the number of Parallel Workers to be used
# print(ncores) # this how many cores are available, and how many you have requested.
# getDoParWorkers()# you can compare with the number of actual workers

# Define parameters
dist <- "log-normal"
params <- list()
censor_params <- list()
censor_dist <- "uniform"
censor_params$start <- 1
censor_params$end <- 8
learner_list = list(base_mars_BinnedIPCW)
params_list <- list(
  list(nprune_values = c(1, 3, 5, 8, 11), degree_values = c(1, 2), useMin = TRUE)
)
time_point = 5
num_obs <- 250
train_prop = 0.7
num_covariates <- 5
params$mean <- 1.5
params$sd <- 0.8
set.seed(2003)
X <- matrix(rnorm(num_obs * num_covariates), nrow = num_obs, ncol = num_covariates)
lc = 2 * sin(pi * X[,1]) + 1.3 * (X[,2]^3) - 1.5 * X[,3] * X[,4] - 0.3 * abs(X[,5]) + 0.6 * X[, 4]
params$lc <- lc
dt <- simulation_data(num_obs, dist, params, censor_dist, censor_params, time_point, X)
Y <- data.frame(E = dt$E, sigma = dt$sigma, observed_time = dt$observed_time)
true_surv = (true_survival_function(dist, time_point, params = params))^exp(lc)
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
df <- combined_df %>%
  group_by(Method) %>%
  summarise(across(1:(ncol(combined_df)-1), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame()


saveRDS(list(results_list = results_list,
             all_layer_testEP = all_layer_testEP,
             combined_results = df), paste0("CM_mars_",dist, "_", num_obs, "_", num_covariates,".rds"))

