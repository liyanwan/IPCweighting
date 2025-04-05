source("~/IPCweighting/dIPCW_ML/sources.R")
source("~/IPCweighting/dIPCW_ML/library.R")

# event time distribution and parameters
dist = "exponential"
params <- list()
params$lambda_base = 0.15

# censoring time distribution and paramters
censor_params <- list()
censor_dist <- "uniform"
censor_params$start <- 1
censor_params$end <- 8

# covariates relevant information
num_obs <- 250
num_covariates <- 25
num_beta <- 25
time_point = 5
train_prop = 0.7
set.seed(2003)
X <- matrix(rnorm(num_obs * num_covariates), nrow = num_obs, ncol = num_covariates)
colnames(X) <- paste0("X", 1:num_covariates)
beta = rnorm(num_beta, mean = 0.1, sd = 0.5)
lc = rowSums(sweep(X,2,beta,"*")) # risk score
true_surv = (true_survival_function(dist, time_point, params = params))^exp(lc) # compute true survival probabilities of all samples
params$lc <- lc # helpful in simulate_data()

# Define the learner list and its parameters
learner_list = list(base_glm_BinnedIPCW)
params_list <- list(
  list(useMin = TRUE)
)

set.seed(69)
# dt is a dataframe consisting of X, event_time, censor_time, observed_time, sigma (delta), E (event status)
dt <- simulation_data(num_obs, dist, params, censor_dist, censor_params, time_point, X)
# M denotes the actual status. 1: event occured; 0: survived more than 5 years.
dt$M = ifelse(dt$event_time<=5, 1, 0)
table(dt$M)
# 0   1 
# 118 132 

# Y is a dataframe which must contain five columns: E, M, event time, sigma, observed_time
Y <- data.frame(E = dt$E, M = dt$M, event_time = dt$event_time,sigma = dt$sigma, observed_time = dt$observed_time)

pp = dIPCW_ML_split(time_point = time_point,  
                    X = X,
                    Y = Y,
                    train_prop = train_prop,  # split X and Y into training and testing data based on train_prop, 
                    measure = "Brier_Score",  # (with stratified (same proportion of E=0,1,NA in both settings)
                    range_intervals = 1,    # Always set to 1 if want Single ML models
                    true_surv = true_surv,   # true survival probabilities of all samples
                    learner_list = learner_list,
                    params_list = params_list,
                    proxy_data = NULL,      # additional data used to fit test_IPCW, usually set it to NULL
                    learner_k = 5,
                    learner_foldid = NULL,  # if learner_foldid = NULL, then use learner_k = 5 to create a 5 folds CV.
                    include_naive = TRUE,   # TRUE if you want to fit naive model
                    include_SA = TRUE,      # TRUE if want survival models
                    include_opt = TRUE,     # TRUE if want binned models
                    max_intervals = 10,     # if include_opt = TRUE, then c(min_intervals, max_intervals) will be the tuning
                    min_intervals = 1).     # grid of bins for binned models

names(pp)
# "results_df"     "all_layer_testEP"    "model_list_with_bin" "model_list_opt"  

pp$results_df
#         Method  C-observed   C-event       AUC     -LogL         BS        OLS Opt_bins
# 1   Single GLM1  0.8893485 0.8545455 0.7926094 2.8388043 0.20547945 0.14890094        1
# 2 OPTBinned GLM  0.8440159 0.8340944 0.9441931 0.3031692 0.10414672 0.02983063        5
# 3     Naive GLM  0.8698482 0.8540707 0.7907240 2.8388043 0.20547945 0.15578624       NA
# 4         coxPH  0.8668588 0.8542308 0.9494721 0.2784630 0.09380957 0.01641162       NA

## C-observed: The C index based on ranking observed times; 
## C-event: The c index based on ranking event times; 
## AUC, -logL, BS: AUC, log-likelihood and Brier score based on actual status (M); OLS: ols error.

pp$all_layer_testEP     ## A data frame of event probabilities for the test data, where each column corresponds to a different method.
#     Single GLM1 OPTBinned GLM    Naive GLM        coxPH
# 1  1.000000e+00  9.999767e-01 1.000000e+00 1.0000000000
# 2  2.220446e-16  1.548420e-03 2.220446e-16 0.0110552223
# 3  2.220446e-16  3.544102e-01 2.220446e-16 0.1671683655
# ......

length(pp$model_list_with_bin) # This has the same length as learner_list
# 1

names(pp$model_list_with_bin[[1]])    ## the information of SINGLE IPCW machine learning model.
# "model"      "bin"    "test_event_pred"   "class_name"

pp$model_list_opt ## The information of the BINNED IPCW machine learning model. It equals to NULL of include_opt = FALSE

names(pp$model_list_with_bin[[1]]) ==  names(pp$model_list_opt[[1]])
# TRUE TRUE TRUE TRUE








