source("~/IPCweighting/dIPCW_ML/sources.R")
source("~/IPCweighting/dIPCW_ML/library.R")

num_obs <- 100
train_prop = 0.7
learner_list = list(base_mars_BinnedIPCW_all)
params_list <- list(
  list(nprune_values = seq(2, 20, by = 2), degree_values = c(1, 2), useMin = TRUE)
)
num_covariates <- 100
num_beta = num_covariates
time_point = 5
rho <- 0.7
cov_matrix <- matrix(0, nrow = num_covariates, ncol = num_covariates)
for (i in 1:num_covariates) {
  for (j in 1:num_covariates) {
    cov_matrix[i, j] = rho^abs(i - j)
  }
}

X = mvrnorm(n = num_obs, mu = rep(0, num_covariates), Sigma = cov_matrix)
colnames(X) <- paste0("X", 1:num_covariates)

multiplier = 1
censor_alpha = 1

beta <- c(rep(0,20), rep(0.1,5),rep(0, 50), rep(-0.1, 5), rep(0, 20))*multiplier#rnorm(num_beta, mean = 0, sd = 0.3)
lc <- rowSums(sweep(X, 2, beta, "*"))
start_iter = 1
end_iter = 40
dist <- "exponential"
params <- list()
params$lambda_base = 0.2
params$lc <- lc
true_surv <- (true_survival_function(dist, time_point, params = params))^exp(lc)
GlobalFunctions = ls(globalenv())

censor_params <- list()
censor_dist <- "Weibull"
censor_params$lambda_base = 0.2
censor_params$alpha = censor_alpha


results_total_list <- foreach(w = start_iter:end_iter, .packages = c("MASS", "dplyr","glmnet", "survival", "caret", "rpart", "earth"),.export=GlobalFunctions) %do% {
  print(w)
  set.seed(20 + w)
  dt <- simulation_data(num_obs, dist, params, censor_dist, censor_params, time_point, X)
  dt$M = ifelse(dt$event_time<=time_point, 1, 0)
  Y <- data.frame(E = dt$E, M = dt$M, event_time = dt$event_time,sigma = dt$sigma, observed_time = dt$observed_time)
  dIPCW_ML_split(time_point = time_point,
                 X = X,
                 Y = Y,
                 train_prop = train_prop,
                 measure = "Brier_Score",
                 range_intervals =1,
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
  summarise(across(1:7, \(x) mean(x, na.rm = TRUE))) %>%  # Compute mean for columns 2 to 7
  as.data.frame()
df
