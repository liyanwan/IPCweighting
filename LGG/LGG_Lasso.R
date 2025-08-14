source("~/IPCweighting/dIPCW_ML/sources.R")
source("~/IPCweighting/dIPCW_ML/library.R")

totalX = readRDS("~/IPCweighting/genes_top300.rds")
totalY = readRDS("~/IPCweighting/Y.rds")
time_point = 5
train_prop = 0.5
learner_list = list(base_glmnet_BinnedIPCW)
params_list <- list(
  list(lambda_values = 10^seq(-1, -5, length = 50), alpha_num = 1, useMin = TRUE)
)

error_occurred <- FALSE
set.seed(0)
totalX = data.frame(totalX)
totalY = data.frame(totalY)
E_factor = factor(totalY$E, levels = c(0, 1, NA), exclude = NULL)
train_index <- createDataPartition(E_factor, p = train_prop, list = FALSE)

results_total_list1 = dIPCW_ML_RD(time_point = time_point, X = totalX[train_index, ], Y = totalY[train_index, ], 
                                  newX = totalX[-train_index, ], newY = totalY[-train_index, ], measure = "Brier_Score", range_intervals = 1, 
                                  true_surv = NULL, learner_list = learner_list, params_list = params_list, proxy_data = NULL, learner_k = 5, learner_foldid = NULL, include_naive = TRUE, include_SA = TRUE, surv_params = list(), include_opt = TRUE, max_intervals = 10, min_intervals = 1)

results_total_list2 = dIPCW_ML_RD(time_point = time_point, X = totalX[-train_index, ], Y = totalY[-train_index, ],
                                  newX = totalX[train_index, ], newY = totalY[train_index, ], measure = "Brier_Score", range_intervals = 1,
                                  true_surv = NULL, learner_list = learner_list, params_list = params_list, proxy_data = NULL, learner_k = 5, learner_foldid = NULL, include_naive = TRUE, include_SA = TRUE, surv_params = list(), include_opt = TRUE, max_intervals = 10, min_intervals = 1)


bind_rows(results_total_list1$results_df, results_total_list2$results_df)%>%
  group_by(Method) %>%
  summarise(across(1:5, \(x) mean(x, na.rm = TRUE))) %>%
  as.data.frame()
