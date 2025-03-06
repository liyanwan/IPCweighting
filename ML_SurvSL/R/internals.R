# Internal functions for SuperLearner package
#
# Created by Eric Polley on 2011-01-03
#

# .SL.require() extends the require() function to add my own error messages
# If a package uses special functions in a formula, may need to use require() instead
.SL.require <- function(package, message = paste('loading required package (', package, ') failed', sep = '')) {
  if(!requireNamespace(package, quietly = FALSE)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

# .createLibrary()
# takes the input from SL.library and sets up the library for SuperLearner()
# SL.library may be a character vector or a list
# return list of prediction algorithms and list of screening algorithms.  If screening algorithms are used, must generate the whichScreen matrix and assign the appropriate row number to the prediction algorithm.  whichScreen[1, ] will always be the full data set, even if no screening algorithms are specified.
# If character vector, only prediction algorithms are allowed.
# for the list, the structure is:
# list(c("predAlg1", "screenAlg1", "screenAlg2", ...), c("predAlg2", "screenAlg1", ...), c("predAlg3"))
# the list contains character vectors, and the scructure of the character vectors is always prediction algorithm first, followed by a list of screening algorithms to be coupled with the prediction algorithm.  If no screening algorithm is given (as in "predAlg3" above) then the algorithm will run on the entire set of variables by default.

.createLibrary <- function(SL.library) {
    if (is.character(SL.library)) {
        k <- length(SL.library)
        whichScreen <- matrix(1, nrow = 1, ncol = k)
        screenAlgorithm <- "All"
        library <- data.frame(predAlgorithm = SL.library, rowScreen = 1, stringsAsFactors=FALSE)
    } else if (is.list(SL.library)) {
        predNames <- sapply(SL.library, FUN = "[", 1)
        NumberScreen <- (sapply(SL.library, FUN = length) - 1)
        if (sum(NumberScreen == 0) > 0) {
            for(ii in which(NumberScreen == 0)) {
                SL.library[[ii]] <- c(SL.library[[ii]], "All")
                NumberScreen[ii] <- 1
            }
        }
        screenAlgorithmFull <- unlist(lapply(SL.library, FUN="[", -1))
        screenAlgorithm <- unique(screenAlgorithmFull)
        
        library <- data.frame(predAlgorithm = rep(predNames, times=NumberScreen), rowScreen = match(screenAlgorithmFull, screenAlgorithm), stringsAsFactors = FALSE)
    } else {
      stop('format for SL.library is not recognized')
    }
    
    out <- list(library = library, screenAlgorithm = screenAlgorithm)
    return(out)
}


# Scans over SL.library and loads required packages.  This step is not necessary since most of the wrapper functions also include require() statements, but this helps detect missing packages early.
# Should this be extended to allow the user to pass additional packages in a character vector to check?  might be helpful for the cluster setting.
# addPackages <- c('MASS', 'class', 'nnet', 'ipred')
# sapply(addPackages, function(x) require(force(x), character.only = TRUE))

.check.SL.library <- function(library, addPackages = NULL) {
    if("dIPCWSL.glmnet" %in% library) .SL.require('glmnet', message = 'You have selected glmnet as a library algorithm but either do not have the glmnet package installed or it can not be loaded')
    if("dIPCWSL.knn" %in% library) .SL.require('yaImpute', message = 'You have selected knn as a library algorithm but either do not have the yaImpute package installed or it can not be loaded')
    if("dIPCWSL.rpart" %in% library) .SL.require('rpart', message = 'You have selected rpart as a library algorithm but either do not have the rpart package installed or it can not be loaded')
    if("dIPCWSL.earth" %in% library) .SL.require('earth', message = 'You have selected earth as a library algorithm but either do not have the earth package installed or it can not be loaded')
    if("dIPCWSL.caret" %in% library) .SL.require('caret', message = 'You have selected caret as a library algorithm but either do not have the caret package installed or it can not be loaded')
#    #removing this check, need to replace with requireNamespace per writing R extensions 1.1.3.1
#    if(!is.null(addPackages)) {
#      sapply(addPackages, function(x) require(force(x), character.only = TRUE))
#    }
    invisible(TRUE)
}
