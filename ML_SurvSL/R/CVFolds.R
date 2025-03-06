#' This function is a modified version of `CVFolds` from the
#' SuperLearner package (Polley et al., 2011).
#' 
# create a list of row numbers for the V-fold cross validation.
# based on sample size N, id, Y, and cvControl
#
# inside cvControl:
# V : number of folds
# stratifyCV : Stratify folds to try and keep proportion constant
# shuffle : should the rows of X, Y be shuffled since the split function is deterministic

CVFolds <- function(N, id, E, cvControl){
  # validRows would be a user specified list of row numbers for the validation sets
    if(!is.null(cvControl$validRows)) {
      return(cvControl$validRows)
    }
    stratifyCV <- cvControl$stratifyCV
    shuffle <- cvControl$shuffle
    V <- cvControl$V
    
    if(!stratifyCV) {
        if(shuffle) {
            if(is.null(id)) {
                validRows <- split(sample(1:N), rep(1:V, length=N))
            } else {
                n.id <- length(unique(id))
                id.split <- split(sample(1:n.id), rep(1:V, length=n.id))
                validRows <- vector("list", V)
                for(v in seq(V)) {
                    validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
                }
            }
        } else {
            if(is.null(id)) {
                validRows <- split(1:N, rep(1:V, length=N))
            } else {
                n.id <- length(unique(id))
                id.split <- split(1:n.id, rep(1:V, length=n.id))
                validRows <- vector("list", V)
                for(v in seq(V)) {
                    validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
                }
            }
        }
    } else {
        if(length(unique(E[!is.na(E)])) != 2) {
            stop("stratifyCV only implemented for binary E (excluding NA values).")
        }
        if(sum(E == 1, na.rm = TRUE) < V | sum(E == 0, na.rm = TRUE) < V) {
            stop("number of (E=1) or (E=0) (excluding NA values) is less than the number of folds.")
        }
        if(shuffle) {
            if(is.null(id)) {
                wiY0 <- which(E == 0)
                wiY1 <- which(E == 1)
                rowsY0 <- split(sample(wiY0), rep(1:V, length=length(wiY0)))
                rowsY1 <- split(sample(wiY1), rep(1:V, length=length(wiY1)))
                validRows <- vector("list", length = V)
                names(validRows) <- paste(seq(V))
                for(vv in seq(V)) {
                   validRows[[vv]] <- c(rowsY0[[vv]], rowsY1[[vv]])
                }
                naRows <- which(is.na(E))
                if(length(naRows) > 0) {
                   naRowsSplit <- split(sample(naRows), rep(1:V, length = length(naRows)))
                   for(vv in seq(V)) {
                       validRows[[vv]] <- c(validRows[[vv]], naRowsSplit[[vv]])
                   }
               }
            } else {
                stop("stratified sampling with id not currently implemented")
            }
        } else {
            if(is.null(id)) {
                within.split <- suppressWarnings(tapply(1:N, INDEX = E, FUN = split, rep(1:V)))
                validRows <- vector("list", length = V)
                names(validRows) <- paste(seq(V))
                for(vv in seq(V)) {
                    validRows[[vv]] <- c(within.split[[1]][[vv]], within.split[[2]][[vv]])
                }
                naRows <- which(is.na(E))
                if(length(naRows) > 0) {
                    naRowsSplit <- split(naRows, rep(1:V, length = length(naRows)))
                    for(vv in seq(V)) {
                        validRows[[vv]] <- c(validRows[[vv]], naRowsSplit[[vv]])
                    }
                }
            } else {
                stop("stratified sampling with id not currently implemented")
            }
        }
    }
    return(validRows)
}
