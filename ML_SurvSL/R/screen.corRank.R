screen.corRank <- function(Y, X, family, method = 'pearson', rank = 15, ...) {
  # if(rank > ncol(X)) {
  #     rank <- ncol(X)
  # } 
  # Don't really need that check, but might want to add a warning message
  E = Y$E
  valid_indices = which(!is.na(E))
  X = X[valid_indices, ]
  E = E[valid_indices]
  listp <- apply(X, 2, function(x, E, method) { 
    ifelse(var(x) <= 0, 1, cor.test(x, y = E, method = method)$p.value)
  }, E = E, method = method)
  whichVariable <- (rank(listp) <= rank)
  return(whichVariable)
}

screen.corRank.l <- function(Y, X, family, method = 'pearson', rank = 20, ...) {
  # if(rank > ncol(X)) {
  #     rank <- ncol(X)
  # } 
  # Don't really need that check, but might want to add a warning message
  E = Y$E
  valid_indices = which(!is.na(E))
  X = X[valid_indices, ]
  E = E[valid_indices]
  listp <- apply(X, 2, function(x, E, method) { 
    ifelse(var(x) <= 0, 1, cor.test(x, y = E, method = method)$p.value)
  }, E = E, method = method)
  whichVariable <- (rank(listp) <= rank)
  return(whichVariable)
}

# screen.corRank <- function(Y, X, family, time_point, method = 'pearson', rank = 2, ...) {
#   # if(rank > ncol(X)) {
#   #     rank <- ncol(X)
#   # } 
#   # Don't really need that check, but might want to add a warning message
#   train_weights = get_IPCW_relevant(data.frame(X, Y), time_point = time_point, return_type = "IPCW")
#   valid_indices = which(train_weights > 0)
#   X = X[valid_indices, ]
#   Y = Y[valid_indices, ]
#   E = Y$E
#   w = train_weights[valid_indices, ]
#   # listp <- apply(X, 2, function(x, Y, method) { 
#   #   ifelse(var(x) <= 0, 1, cor.test(x, y = Y, method = method)$p.value)
#   # }, Y = Y, method = method)
#   listp <- apply(X, 2, function(x, E, w, method) { 
#     if (var(x) <= 0) {1} else {
#       weighted.cor.test(x, E = E, w = w, method = method)$p.value
#     }
#   }, E = E, w = w, method = method)
#   whichVariable <- (rank(listp) <= rank)
#   return(whichVariable)
# }