# control functions for dIPCWSuperLearner()
#
# Based on the SuperLearner package (Polley et al., 2011) - URL: https://cran.r-project.org/web/packages/SuperLearner/
#
dIPCWSuperLearner.control <- function(saveFitLibrary = TRUE,
                                 saveCVFitLibrary = FALSE,
                                 trimLogit = 0.001) {
  if(trimLogit > 0.5) {
    warning('trimLogit must be less than 0.5, will replace with trimLogit = 0.001')
    trimLogit <- 0.001
  }
  list(saveFitLibrary = saveFitLibrary, trimLogit = trimLogit, saveCVFitLibrary = saveCVFitLibrary)
}

dIPCWSuperLearner.CV.control <- function(V = 10L, stratifyCV = FALSE, shuffle = TRUE, validRows = NULL){
  # make sure V is an integer
  V <- as.integer(V)
  if(!is.null(validRows)) {
    if(!is.list(validRows)) {
      stop('validRows must be a list of length V containing the row numbers for the corresponding validation set')
    }
    if(!identical(V, length(validRows))) {
      stop('V and length(validRows) must be identical')
    }
  }
  list(V = V, stratifyCV = stratifyCV, shuffle = shuffle, validRows = validRows)
}
