source("~/IPCweighting/dIPCW_ML/sources.R")
source("~/IPCweighting/dIPCW_ML/library.R")

# ncores = detectCores()
# registerDoParallel(cores=ncores-3)# Shows the number of Parallel Workers to be used
# print(ncores) # this how many cores are available, and how many you have requested.
# getDoParWorkers()# you can compare with the number of actual workers

# Define parameters
dist <- "exponential"
params <- list()
params$lambda_base = 0.15
censor_params <- list()
censor_dist <- "uniform"
censor_params$start <- 1
censor_params$end <- 8

num_obs <- 250
train_prop = 0.7
learner_list = list(base_glmnet_BinnedIPCW)
params_list <- list(
  list(lambda_values = 10^seq(-1, -5, length = 50), alpha_num = 1, useMin = TRUE)
)
num_covariates <- 100
num_beta = 25

set.seed(2003)
time_point = 5
X <- matrix(rnorm(num_obs * num_covariates), nrow = num_obs, ncol = num_covariates)
colnames(X) <- paste0("X", 1:num_covariates)
start_iter = 1
end_iter = 2
# Used to check if there will be infinite event time in 1000 simulations. (applied to log-normal only)
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

results_total_list <- foreach(w = start_iter:end_iter, 
                              .packages = c("MASS", "dplyr","glmnet", "survival", "caret", "rpart", "earth"),
                              .export=GlobalFunctions) %do% {
    set.seed(20 + w)
    beta <- c(rnorm(num_beta, mean = 0, sd = 0.3), rep(0, num_covariates-num_beta))
    lc <- rowSums(sweep(X, 2, beta, "*"))
    params$lc <- lc
    dt <- simulation_data(num_obs, dist, params, censor_dist, censor_params, time_point, X)
    dt$M = ifelse(dt$event_time<=5, 1, 0)
    table(dt$M)
    Y <- data.frame(E = dt$E, M = dt$M, event_time = dt$event_time,sigma = dt$sigma, observed_time = dt$observed_time)
    true_surv <- (true_survival_function(dist, time_point, params = params))^exp(lc)
    set.seed(20+w)
    dIPCW_ML_split(time_point = time_point,
                   X = X,
                   Y = Y,
                   train_prop = train_prop,
                   measure = "Brier_Score",
                   range_intervals = 1,
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

saveRDS(list(all_layer_testEP = all_layer_testEP,
             results_list = results_list,
             combined_results = df),
        paste0("BinnedHDall_", dist, "_", num_obs, "_", num_covariates, ".rds"))
