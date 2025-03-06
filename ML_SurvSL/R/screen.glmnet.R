screen.glmnet <- function(Y, X, family, alpha = 1, minscreen = 2, nfolds = 10, lambda = 10^seq(-1, -5, length = 50), threshold = NULL,  ...) {
  .SL.require('glmnet')
  E = Y$E
  valid_indices = which(!is.na(E))
  X = X[valid_indices, ]
  E = E[valid_indices]
  if(!is.matrix(X)) {
    X <- model.matrix(~ -1 + ., X)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = E, type.measure = 'deviance', nfolds = nfolds, family = family$family, alpha = alpha, lambda = lambda)
  whichVariable <- (as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1] != 0)
  # the [-1] removes the intercept
  if (sum(whichVariable) < minscreen) {
    warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
    sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta), 2, function(x) sum((x != 0)))
    newCut <- which.max(sumCoef >= minscreen)
    whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[, newCut] != 0)
  }
  if((!is.null(threshold)) && (sum(whichVariable) > threshold)){
    coef_values <- as.numeric(coef(fitCV$glmnet.fit, s = fitCV$lambda.min))[-1]
    abs_coef <- abs(coef_values)
    whichVariable <- (rank(-abs_coef) <= threshold)
  }
  return(whichVariable)
}