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
censor_params <- list()
censor_dist <- "uniform"
censor_params$start <- 1
censor_params$end <- 8
learner_list = list(base_CTree_BinnedIPCW)
params_list <- list(
  list(ccp_alpha_values = 2 * 10^seq(-1, -4, length = 40), useMin = TRUE)
)

num_obs <- 250
time_point = 5
train_prop = 0.7
num_covariates <- 15
num_beta <- 15
rho <- 0.8
cov_matrix <- matrix(0, nrow = num_covariates, ncol = num_covariates)
for (i in 1:num_covariates) {
  for (j in 1:num_covariates) {
    cov_matrix[i, j] = rho^abs(i - j)
  }
}

set.seed(2003)
X = mvrnorm(n = num_obs, mu = rep(0, num_covariates), Sigma = cov_matrix)
X_model = mvrnorm(n = num_obs, mu = rep(0, num_covariates), Sigma = cov_matrix)
colnames(X_model) = paste0("X", 1:num_covariates)
colnames(X) = paste0("X", 1:num_covariates)
beta = rnorm(num_beta, mean = 0.6, sd = 0.01)
Y_continuous = generate_Y_continuous(X_model, c(1:num_covariates), beta)+rnorm(num_obs)
data_tree = as.data.frame(X_model)
data_tree$Y_continuous = Y_continuous
tree_model_fX = rpart(Y_continuous ~ ., data = data_tree,  control = rpart.control(cp = 0.03))
f_X = predict(tree_model_fX, newdata =as.data.frame(X))
lc = f_X
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
             combined_results = df), paste0("CM_trees_",dist, "_", num_obs, "_", num_covariates,".rds"))

