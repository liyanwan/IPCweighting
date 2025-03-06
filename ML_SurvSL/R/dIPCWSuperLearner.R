#  SuperLearner for right-censored data with discrete-time interval
#
dIPCWSuperLearner <- function(Y, X, time_point, newX = NULL, family = binomial(), SL.library, measure,
                        method = 'method.glmnet', meta_learner_params = NULL, id = NULL, verbose = FALSE,
                        control = list(), cvControl = list(), env = parent.frame()) {
   # Y is a dataframe with column observed_time, E, and sigma.
   # Begin timing how long SuperLearner takes to execute.
   time_start = proc.time()
   if(is.matrix(X)){
     X = as.data.frame(X)
     Y = as.data.frame(Y)
     newX = as.data.frame(newX)
   }

   if (is.character(method)) {
       if (exists(method, mode = 'list')) {
           method <- get(method, mode = 'list')
       } else if (exists(method, mode = 'function')) {
           method <- get(method, mode = 'function')()
       }
   } else if (is.function(method)) {
       method <- method()
   }
   if(!is.list(method)) {
       stop("method is not in the appropriate format. Check out help('method.template')")
   }
   if(!is.null(method$require)) {
       sapply(method$require, function(x) require(force(x), character.only = TRUE))
   }
   # get defaults for controls and make sure in correct format
   control <- do.call('dIPCWSuperLearner.control', control)
   cvControl <- do.call('dIPCWSuperLearner.CV.control', cvControl)

   library <- .createLibrary(SL.library)
   .check.SL.library(library = c(unique(library$library$predAlgorithm), library$screenAlgorithm))

   call <- match.call(expand.dots = TRUE)
   # should we be checking X and newX for data.frame?
   # data.frame not required, but most of the built-in wrappers assume a data.frame
   if(!inherits(X, 'data.frame')) message('X is not a data frame. Check the algorithms in SL.library to make sure they are compatible with non data.frame inputs')
   varNames <- colnames(X)
   N <- dim(X)[1L]
   p <- dim(X)[2L]
   k <- nrow(library$library)
   kScreen <- length(library$screenAlgorithm)
   
   if(p < 2 & !identical(library$screenAlgorithm, "All")) {
       warning('Screening algorithms specified in combination with single-column X.')
   }

   # put fitLibrary in it's own environment to locate later
   fitLibEnv <- new.env()
   

   # if newX is missing, use X
   if(is.null(newX)) {
       newX <- X
   }
   # Are these checks still required?
   if(!identical(colnames(X), colnames(newX))) {
       stop("The variable names and order in newX must be identical to the variable names and order in X")
   }
   if (sum(is.na(X)) > 0 | sum(is.na(newX)) > 0) {
       stop("missing data is currently not supported. Check X, and newX for missing values")
   }
   
   # family can be either character or function, so these lines put everything together (code from glm())
   if(family$family != "binomial"){
       stop("Gaussian is currently not supported.")
   }
   if(is.character(family))
       family <- get(family, mode="function", envir=parent.frame())
   if(is.function(family))
       family <- family()
   if (is.null(family$family)) {
       print(family)
       stop("'family' not recognized")
   }
   # create CV folds, CVFolds is able to handle NA values in E. default id = NULL
   E = Y$E
   validRows <- CVFolds(N = N, id = id, E = E, cvControl = cvControl)

   # test id
   if(is.null(id)) {
       id <- seq(N)
   }
   if(!identical(length(id), N)) {
       stop("id vector must have the same dimension as E")
   }
   
   time_train_start = proc.time()
   
   # cross-validation step
   cvout = dIPCWSL.crossValStep(validRows = validRows,
                            Y = Y, dataX = X, time_point = time_point, 
                            id = id, family = family, library = library, kScreen = kScreen, k = k,
                            p = p, env = env, saveCVFitLibrary = control$saveCVFitLibrary)
   zz = cvout$Z
   if(control$saveCVFitLibrary){
     cvFitLibrary <- cvout$model_out
   }else{
     cvFitLibrary <- NULL
   }
   libraryNames <- paste(library$library$predAlgorithm, library$screenAlgorithm[library$library$rowScreen], sep="_")

   q = length(libraryNames)
   assign('fitLibrary', vector('list', length = q), envir = fitLibEnv)
   assign('libraryNames', libraryNames, envir = fitLibEnv)
   evalq(names(fitLibrary) <- libraryNames, envir = fitLibEnv)
   
   # errors* records if an algorithm stops either in the CV step and/or in full data
   errorsInCVLibrary <- rep(0, q)
   errorsInLibrary <- rep(0, q)
   
   # Check for errors. If any algorithms had errors, replace entire column with
   # 0 even if error is only in one fold.
   errorsInCVLibrary <- apply(zz, 2, function(x) anyNA(x))
   if (sum(errorsInCVLibrary) > 0) {
     zz[, as.logical(errorsInCVLibrary)] <- 0
   }
   if (all(zz == 0)) {
       stop("All algorithms dropped from library")
   }
   foldid = rep(NA, nrow(Y))
   for (i in seq_along(validRows)) {
     foldid[validRows[[i]]] = i
   }
   time_train = proc.time() - time_train_start

   # now fit all algorithms in library on entire learning data set and predict on newX
   m <- dim(newX)[1L]
   predE <- matrix(NA, nrow = m, ncol = q)
   # whichScreen <- matrix(NA, nrow = kScreen, ncol = p)

   .screenFun <- function(fun, list) {
       screen_fn = get(fun, envir = env)
       testScreen <- try(do.call(screen_fn, list))
       if (inherits(testScreen, "try-error")) {
           warning(paste("replacing failed screening algorithm,", fun, ", with All() in full data", "\n "))
           out <- rep(TRUE, ncol(list$X))
       } else {
           out <- testScreen
       }
       return(out)
   }

   time_predict_start = proc.time()

   whichScreen <- sapply(library$screenAlgorithm, FUN = .screenFun, list = list(Y = Y, X = X, family = family, id = id), simplify = FALSE)
   whichScreen <- do.call(rbind, whichScreen)
   
   # Fit each candidate in library to X
   .predFun <- function(index, lib, Y, dataX, newX, time_point, whichScreen, family, id, verbose, control, libraryNames) {
     pred_fn = get(lib$predAlgorithm[index], envir = env)
     testAlg <- try(do.call(pred_fn, list(Y = Y,
                                          X = subset(dataX,
                                                     select = whichScreen[lib$rowScreen[index], ], drop=FALSE),
                                          newX = subset(newX, select = whichScreen[lib$rowScreen[index], ], drop=FALSE),
                                          family = family, id = id, time_point = time_point)))
     if (inherits(testAlg, "try-error")) {
       warning(paste("Error in algorithm", lib$predAlgorithm[index], " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n" ))
       out <- rep.int(NA, times = nrow(newX))
     } else {
       out <- testAlg$pred
       if (control$saveFitLibrary) {
         eval(bquote(fitLibrary[[.(index)]] <- .(testAlg$fit)), envir = fitLibEnv)
       }
     }
     if (verbose) {
       message(paste("full", libraryNames[index]))
     }
     invisible(out)
   }
   
   predE <- do.call('cbind', lapply(seq(k), FUN = .predFun,
                                    lib = library$library, Y = Y, dataX = X,
                                    newX = newX, time_point = time_point, whichScreen = whichScreen,
                                    family = family, id = id,
                                    verbose = verbose, control = control,
                                    libraryNames = libraryNames))
   
   # assign('fitLibrary', foo$fitLibrary, envir = fitLibEnv)
   time_predict = proc.time() - time_predict_start
   
   # check for errors
   errorsInLibrary <- apply(predE, 2, function(algorithm) anyNA(algorithm))
   if (sum(errorsInLibrary) > 0) {
       if (sum(coef[as.logical(errorsInLibrary)]) > 0) {
           warning(paste0("Re-running estimation of coefficients removing failed algorithm(s)\n",
                          "Original coefficients are: \n"))
           zz[, as.logical(errorsInLibrary)] <- 0
           if (all(zz == 0)) {
               stop("All algorithms dropped from library")
           }
       } else {
           warning("Coefficients already 0 for all failed algorithm(s)")
       }
   }
   
   # Compute super learner predictions on newX.
   trainobsWeights = get_IPCW_relevant(observed_time = Y$observed_time, sigma = Y$sigma, 
                                       time_point = time_point, return_type = "IPCW")
   getCoef_bin <- method$computeCoef(Z = zz, Y = Y, libraryNames = libraryNames, trainobsWeights = trainobsWeights, 
                                     folds = foldid, time_point = time_point, measure = measure,
                                     control = control, verbose = verbose, alpha_num = meta_learner_params$alpha_num,
                                     lambda_grid = meta_learner_params$lambda_values, useMin = meta_learner_params$useMin,
                                     errorsInLibrary = errorsInCVLibrary)
   coef = getCoef_bin$coef
   if(method$require == "glmnet"){
     intercept = getCoef_bin$intercept
     coef_total = c(intercept, coef)
   } else {
     coef_total = coef
   }
   getPred <- method$computePred(predE = predE, coef = coef_total, control = control)

   # Add names of algorithms to the predictions.
   colnames(predE) <- libraryNames

   # Clean up when errors in library.
   if(sum(errorsInCVLibrary) > 0) {
       getCoef$cvRisk[as.logical(errorsInCVLibrary)] <- NA
   }

   # Finish timing the full SuperLearner execution.
   time_end = proc.time()

   # Compile execution times.
   times = list(everything = time_end - time_start,
                train = time_train,
                predict = time_predict)

   # Put everything together in a list.
   out <- list(
       call = call,
       libraryNames = libraryNames,
       cvFitLibrary = cvFitLibrary,
       SL.library = library,
       SL.predict = getPred,
       coef = coef_total,
       library.predict = predE,
       Z = zz,
       cvRisk = getCoef_bin$cvRisk,
       MLCVRisk = getCoef_bin$ML_CVRisk,
       family = family,
       fitLibrary = get('fitLibrary', envir = fitLibEnv),
       varNames = varNames,
       validRows = validRows,
       method = method,
       whichScreen = whichScreen,
       control = control,
       cvControl = cvControl,
       errorsInCVLibrary = errorsInCVLibrary,
       errorsInLibrary = errorsInLibrary,
       env = env,
       times = times
   )
   class(out) <- c("dIPCWSuperLearner")
   return(out)
}
