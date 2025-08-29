# Binned Inverse Probability of Censoring Weighting

## Overview
This repository contains the implementation and resources for the project **"Integrating Binning with Inverse Probability of Censoring Weighting for Improved Risk Prediction with Machine Learning"**. The goal of this project is to introduce the **BIPCW** technique and evaluate its performance against other approaches when applied to various classification machine learning models under Cox proportional hazards (CoxPH) settings. The repository includes simulation frameworks that support multiple baseline survival functions and censoring time distributions.

## Example

```r
BIPCW_ML(
    time_point,
    X,
    Y,
    newX = NULL,
    newY = NULL,
    measure = c('C_index', 'Log_Likelihood_Neg', "Brier_Score")
    range_intervals,
    true_surv = NULL,
    learner_list,
    params_list,
    proxy_data = NULL,
    learner_k = 5,
    learner_foldid = NULL,
    include_naive = FALSE,
    include_SA = FALSE,
    surv_params = list(),
    include_opt = FALSE,
    max_intervals = 10,
    min_intervals = 1
)
```

### Arguments

- time_point: Time horizon.  
- X: Covariates matrix for training.  
- Y: Data frame with columns: `E`, `M`, `event time`, `sigma`, `observed_time` (training data).  
- newX, newY: Test data.  
- measure: Tuning selection metric. One of c('C_index', 'Log_Likelihood_Neg', 'Brier_Score')  
- range_intervals: Grid of positive integers for tuning `bins`.  
- true_surv: True survival probabilities for `newX`.  
- learner_list: List of learner functions.  
- params_list: Nested list of parameters for each learner.  
- proxy_data: Additional data for computing `test_IPCW` (set to `NULL` in experiments).  
- learner_k, learner_foldid: Cross-validation settings. If `learner_foldid` is `NULL`, folds are generated using `learner_k`. All models share the same folds.  
- include_naive: `TRUE` to fit naive models.  
- include_SA: `TRUE` to fit survival models corresponding to learners.  
- surv_params: Default `list()`; only set when fitting glm, glmnet, MARS and Neural Network. The length of surv_params should be equal to the length of params_list.
- include_opt: `TRUE` to include binned IPCW models. Adjust bin range with `max_intervals` and `min_intervals`.



## Available Distributions

### Baseline Survival Functions
- Exponential, Weibull, Log-normal, Log-logistic  

### Censor Time Distributions
- Exponential, Weibull, Log-normal, Log-logistic, Uniform, Gamma  

## Simulation Experiments
We adapt several machine learning models to the BIPCW framework, including GLM, GLMnet, Classification Tree, MARS, and Neural Networks.

The following scripts implement the experiments:
 - `Sim_Experiments_Lasso.R`
 - `Sim_Experiments_Tree.R`
 - `Sim_Experiments_MARS.R`

These experiments compare the performance of: Traditional survival analysis models, BIPCW + ML, IPCW + ML and Naïve ML

Performance is evaluated using the concordance index (C-index), AUC, negative log-likelihood, and Brier score, based on predicted 5-year overall survival probabilities.

Modify **`multiplier`** and **`alpha`** in scripts to change data generation settings.


## TCGA Data Application 
 - `LGG_ExtractInfo.R`: Extracts and preprocesses gene expression profiles from TCGA.  The IDH subtype for each patient is provided in `IDH.csv`, which was extracted by querying TCGA using `lgg.gbm.subtype <- TCGAquery_subtype(tumor = "lgg")`.
 - `LGG_Lasso.R`, `LGG_Tree.R`, `LGG_MARS.R`, `LGG_nn.R`: Each script compares three model variants (except Naive) to compare their performance for TCGA LGG data.
