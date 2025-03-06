dIPCWSL.crossValStep <- function(validRows, Y, dataX, time_point, id, family, library, kScreen, k, p, env, saveCVFitLibrary) {
  num_folds = length(validRows)
  if(saveCVFitLibrary){
    model_out <- vector(mode = "list", length = k)
  } else {
    model_out <- NULL
  }
  results = lapply(seq(k), function(s){
    pred_fn <- get(library$library$predAlgorithm[s], envir = env)
    candidate_col = matrix(NA, nrow = nrow(dataX), ncol = 1)
    for(i in 1:num_folds){
      validation_idx = validRows[[i]]
      X_train = dataX[-validation_idx, ]
      y_train = Y[-validation_idx, ]
      X_validation = dataX[validation_idx, ]
      tempId = id[-validation_idx]
      tempWhichScreen <- matrix(NA, nrow = kScreen, ncol = p)
      for(e in seq(kScreen)) {
        screen_fn = get(library$screenAlgorithm[e], envir = env)
        testScreen <- try(do.call(screen_fn, list(Y = y_train, X = X_train, family = family, id = tempId, time_point = time_point)))
        if(inherits(testScreen, "try-error")) {
          warning(paste("replacing failed screening algorithm,", library$screenAlgorithm[e], ", with All()", "\n "))
          tempWhichScreen[e, ] <- TRUE
        } else {
          tempWhichScreen[e, ] <- testScreen
        }
      }
      testAlg <- try(do.call(pred_fn, list(
        Y = y_train, 
        X = subset(X_train, select = tempWhichScreen[library$library$rowScreen[s], ], drop = FALSE), 
        newX = subset(X_validation, select = tempWhichScreen[library$library$rowScreen[s], ], drop = FALSE), 
        family = family, 
        id = tempId, 
        time_point = time_point
      )))
      if (inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", library$library$predAlgorithm[s], 
                      "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      } else {
        preds = as.matrix(testAlg$pred)
        candidate_col[validation_idx,] = preds
        if(saveCVFitLibrary){
          model_out[[s]][[i]] = testAlg$fit
          }
        }
      }
    candidate_col
    })
  Z = do.call(cbind, results)
  return(list(Z = Z, model_out = model_out))
  }