source("~/IPCweighting/dIPCW_ML/sources.R")
source("~/IPCweighting/dIPCW_ML/library.R")

# ncores = 25
# registerDoParallel(cores=ncores)# Shows the number of Parallel Workers to be used
# print(ncores) # this how many cores are available, and how many you have requested.
# getDoParWorkers()# you can compare with the number of actual workers

# Define parameters
dist = "exponential"
params <- list()
censor_params <- list()
censor_dist <- "uniform"
censor_params$start <- 1
censor_params$end <- 8
learner_list = list(base_glm_BinnedIPCW)
params_list <- list(
  list(useMin = TRUE)
)
num_obs <- 300
num_covariates <- 15
num_beta <- 15
#params$mean <- 1.5
#params$sd <- 0.8
params$lambda_base = 0.15
time_point = 5
train_prop = 0.7
set.seed(2003)
X <- matrix(rnorm(num_obs * num_covariates), nrow = num_obs, ncol = num_covariates)
colnames(X) <- paste0("X", 1:num_covariates)
beta = rnorm(num_beta, mean = 0.65, sd = 0.01)
lc = rowSums(sweep(X,2,beta,"*"))
true_surv = (true_survival_function(dist, time_point, params = params))^exp(lc)
params$lc <- lc
dt <- simulation_data(num_obs, dist, params, censor_dist, censor_params, time_point, X)
Y <- data.frame(E = dt$E, sigma = dt$sigma, observed_time = dt$observed_time)

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
                                               range_intervals = 8,
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
             combined_results = df), paste0("CM_ph_",dist, "_", num_obs, "_", num_covariates,".rds"))

