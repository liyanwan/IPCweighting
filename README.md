# Inverse Probability of Censoring Weighting

## Experiments
The `Comparison_models.Rmd` file simulates a right-censored dataset and compares different estimation methods. The following arguments can be specified for the simulations:
 `dist` `censor_dist`: Distribution of event_time and censor_time.
 `seed`: Number of seeds used in the simulations.
 `num_obs`: Number of observations in each dataset.
 `time_point`: A predefined time point. (Default: 5 years)
We assume the true survival function follows the model: S(t|X_i) = (S_0(t))^exp(X * beta)
 The file prints the OLS_Error(to be continued) between the true survival function and the estimated survival function, corresponding to the following models:
 1. Naive estimate
 2. IPCW (1 bin) estimate
 3. True baseline survival function and estimated beta (from COX Proportional Hazard model)
 4. Survival function with estimated baseline survival function and estimated beta (from COX Proportional Hazard model)

 Available Distribution for event_time:
   - Exponential
   - Weibull
   - Log-logistic
   - Log-Normal
Available Distribution for censor_time:
   - Exponential
   - Weibull
   - Log-logistic
   - Log-normal
   - Uniform

The `binned_experiments.Rmd` file simulates a right-censored dataset, finds the optimal number of bins, and compares the performance of binned IPCW with one dataset and that of binned IPCW with n datasets (`n` means the number of bins). The following arguments can be specified:
   `dist` `censor_dist`: Distribution of event_time and censor_time.
   `seed`: Number of seeds used in the simulations.
   `num_obs`: Number of observations in each dataset.
   `time_point`: A predefined time point. (Default: 5 years)
   `max_intervals`: Maximum number of intervals. (Default: 80)
The file generates two plots: 
    1. OLS_Error(to be continued) vs. the number of intervals for one dataset case.
    2. OLS_prints vs. the number of intervals for `n` datasets case.

 

