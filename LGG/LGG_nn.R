totalX = readRDS("~/IPCweighting/genes_top300.rds")
totalY = readRDS("~/IPCweighting/Y100.rds")
time_point = 7.5
train_prop = 0.5

set.seed(0)
model = dNNmodel(
  units = c(64, 32,1),
  activation = c("elu", "elu","sigmoid"),
  input_shape = 300
)

learner_list = list(base_deepnn_BinnedIPCW)
params_list <- list(
  list(alpha = c(0.1, 0.4, 0.7, 0.9), lambda = c(1e-3, 1e-2, 0.1, 0.3, 1, 3, 10), lr_rate = c(0.01, 0.001, 0.0001),
       model = model,
       useMin = TRUE,
       batch_size = 64)
)

error_occurred <- FALSE
totalY$E = ifelse(totalY$observed_time <= time_point & totalY$sigma == 1, 1, 
                  ifelse(totalY$observed_time > time_point, 0, NA))
set.seed(0)
totalX = data.frame(totalX)
totalY = data.frame(totalY)
E_factor = factor(totalY$E, levels = c(0, 1, NA), exclude = NULL)
train_index <- createDataPartition(E_factor, p = train_prop, list = FALSE)
GlobalFunctions = ls(globalenv())

trainE_stratified = ifelse(is.na(totalY[train_index, ]$E), "missing", as.character(totalY[train_index, ]$E))
train_E_factor = factor(trainE_stratified)
cv_folds <- createFolds(train_E_factor, k = 5, returnTrain = FALSE)
learner_foldid_1 = rep(NA, nrow(totalY[train_index, ]))
for (i in seq_along(cv_folds)) {
  learner_foldid_1[cv_folds[[i]]] = i
}

trainE_stratified = ifelse(is.na(totalY[-train_index, ]$E), "missing", as.character(totalY[-train_index, ]$E))
train_E_factor = factor(trainE_stratified)
cv_folds <- createFolds(train_E_factor, k = 5, returnTrain = FALSE)
learner_foldid_2 = rep(NA, nrow(totalY[-train_index, ]))
for (i in seq_along(cv_folds)) {
  learner_foldid_2[cv_folds[[i]]] = i
}




results_total_list1_time7.5 = lapply(c(1, 2, 5,8), function(bin){
  print(bin)
  if(bin == 1){
    dIPCW_ML_RD(time_point = time_point, X = totalX[train_index, ], Y = totalY[train_index, ],
                newX = totalX[-train_index, ], newY = totalY[-train_index, ], measure = "Brier_Score",
                range_intervals = bin, true_surv = NULL, learner_list = learner_list, params_list = params_list,
                proxy_data = NULL, learner_k = 5, learner_foldid = learner_foldid_1, include_naive = FALSE,
                include_SA = TRUE, surv_params = list(), include_opt = TRUE, max_intervals = 10, min_intervals = 1)
  } else {
    dIPCW_ML_RD(time_point = time_point, X = totalX[train_index, ], Y = totalY[train_index, ],
                newX = totalX[-train_index, ], newY = totalY[-train_index, ],
                measure = "Brier_Score", range_intervals = bin, true_surv = NULL,
                learner_list = learner_list, params_list = params_list, proxy_data = NULL,
                learner_k = 5, learner_foldid = learner_foldid_1, include_naive = FALSE, include_SA = FALSE,
                surv_params = list(), include_opt = FALSE, max_intervals = 10, min_intervals = 1)
  }
})



results_total_list2_time7.5 = lapply(c(1,2,5,8), function(bin){
  if(bin == 1){
    print(bin)
    dIPCW_ML_RD(time_point = time_point, X = totalX[-train_index, ], Y = totalY[-train_index, ],
                newX = totalX[train_index, ], newY = totalY[train_index, ], measure = "Brier_Score", range_intervals = bin,
                true_surv = NULL, learner_list = learner_list, params_list = params_list,
                proxy_data = NULL, learner_k = 5, learner_foldid = learner_foldid_2, include_naive = FALSE,
                include_SA = TRUE, surv_params = list(), include_opt = TRUE, max_intervals = 10, min_intervals = 1)
  } else {
    print(bin)
    dIPCW_ML_RD(time_point = time_point, X = totalX[-train_index, ], Y = totalY[-train_index, ],
                newX = totalX[train_index, ], newY = totalY[train_index, ], measure = "Brier_Score", range_intervals = bin,
                true_surv = NULL, learner_list = learner_list, params_list = params_list,
                proxy_data = NULL, learner_k = 5, learner_foldid = learner_foldid_2, include_naive = FALSE,
                include_SA = FALSE, surv_params = list(), include_opt = FALSE, max_intervals = 10, min_intervals = 1)
  }
})


