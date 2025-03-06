screen.corP <- function(Y, X, family, obsWeights, id, method = 'pearson', minPvalue = 0.1, minscreen = 2, ...) {
  E = Y$E
  valid_indices = which(!is.na(E))
  X = X[valid_indices, ]
  E = E[valid_indices]
  listp <- apply(X, 2, function(x, E, method) { 
    ifelse(var(x) <= 0, 1, cor.test(x, y = E, method = method)$p.value)
  }, E = E, method = method)
  whichVariable <- (listp <= minPvalue)
  which(whichVariable)
  if (sum(whichVariable) < minscreen) {
    warning('number of variables with p value less than minPvalue is less than minscreen')
    whichVariable[rank(listp) <= minscreen] <- TRUE
  }
  return(whichVariable)
}

screen.corP.s <- function(Y, X, family, obsWeights, id, method = 'pearson', minPvalue = 0.02, minscreen = 2, ...) {
  E = Y$E
  valid_indices = which(!is.na(E))
  X = X[valid_indices, ]
  E = E[valid_indices]
  listp <- apply(X, 2, function(x, E, method) { 
    ifelse(var(x) <= 0, 1, cor.test(x, y = E, method = method)$p.value)
  }, E = E, method = method)
  whichVariable <- (listp <= minPvalue)
  which(whichVariable)
  if (sum(whichVariable) < minscreen) {
    warning('number of variables with p value less than minPvalue is less than minscreen')
    whichVariable[rank(listp) <= minscreen] <- TRUE
  }
  return(whichVariable)
}


# screen.corP <- function(E, train_X, family, time_point, id, method = 'pearson', minPvalue = 0.1, minscreen = 2, ...) {
#   train_weights = get_IPCW_relevant(data.frame(train_X, train_Y), time_point = time_point, return_type = "IPCW")
#   valid_indices = which(train_weights > 0)
#   train_X = train_X[valid_indices, ]
#   E = E[valid_indices]
#   # Y = Y[valid_indices, ]
#   # E = Y$E
#   w = train_weights[valid_indices]
#   listp_weighted <- apply(train_X, 2, function(x, E, w, method) { 
#     if (var(x) <= 0) {1} else {
#       weighted.cor.test(x, E = E, w = w, method = method, alternative = "two.sided")$p.value
#     }
#   }, E = E, w = w, method = method)
#   whichVariable_weighted <- (listp <= minPvalue)
#   which(whichVariable_weighted)
#   
#   if (sum(whichVariable) < minscreen) {
#     warning('number of variables with p value less than minPvalue is less than minscreen')
#     whichVariable[rank(listp) <= minscreen] <- TRUE
#   }
#   return(whichVariable)
# }

# weighted.cor.test <- function(x, E, w, alternative = c("two.sided", "less", "greater"),
#                               method = c("pearson", "kendall", "spearman"), exact = NULL,
#                               conf.level = 0.95, continuity = FALSE, ...) { #NO NA in E
#   alternative <- match.arg(alternative)
#   # method <- match.arg(method)
#   DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(E)))
#   if(is.null(w)) {
#     w <- rep(1, length(x))
#   } else {
#     if(length(w) != length(x)) {
#       stop("'w' must be the same length as 'x' and 'E'.")
#     }
#   }
#   if(length(x) != length(E))
#     stop("'x' and 'E' must have the same length")
#   if(!is.numeric(x)) stop("'x' must be a numeric vector")
#   if(!is.numeric(E)) stop("'E' must be a numeric vector")
#   OK <- complete.cases(x, E)
#   x <- x[OK]
#   E <- E[OK]
#   w <- w[OK]
#   n <- length(x)
#   
#   NVAL <- 0
#   conf.int <- FALSE
#   
#   if(method == "pearson") {
#     if(n < 3L)
#       stop("not enough finite observations")
#     method <- "Pearson's product-moment correlation"
#     names(NVAL) <- "correlation"
#     w_sum <- sum(w)
#     r = weights::wtd.cor(x = x, y = E, weight = w)
#     # r <- cor(x, E)
#     fit <- lm(E ~ x, weights = w)
#     smry <- summary(fit)
#     coefs <- smry$coefficients
#     slope_est <- coefs[2, 1]  # estimate of slope
#     slope_se  <- coefs[2, 2]  # std error
#     t_value   <- coefs[2, 3]  # t-statistic
#     p_value   <- coefs[2, 4]  # p-value (two-sided)
#     df_resid <- smry$df[2]
#     ESTIMATE <- c(cor = r)
#     PARAMETER <- c(df = df_resid)
#     STATISTIC <- c(t = t_value) # c(t = sqrt(df) * r / sqrt(1 - r^2))
#     PVAL <- switch(alternative,
#                    "less" = pt(STATISTIC, df_resid),
#                    "greater" = pt(STATISTIC, df_resid, lower.tail=FALSE),
#                    "two.sided" = p_value)
#                    # "two.sided" = 2 * min(pt(STATISTIC, df_resid),
#                    #                       pt(STATISTIC, df_resid, lower.tail=FALSE)))
#     
#     if(n > 3) { 
#       if(!missing(conf.level) &&
#          (length(conf.level) != 1 || !is.finite(conf.level) ||
#           conf.level < 0 || conf.level > 1))
#         stop("'conf.level' must be a single number between 0 and 1")
#       conf.int <- TRUE
#       z <- atanh(r[1])
#       sigma <- 1 / sqrt(w_sum - 3)
#       cint <-
#         switch(alternative,
#                less = c(-Inf, z + sigma * qnorm(conf.level)),
#                greater = c(z - sigma * qnorm(conf.level), Inf),
#                two.sided = z + c(-1, 1) * sigma * qnorm((1 + conf.level) / 2))
#       cint <- tanh(cint)
#       attr(cint, "conf.level") <- conf.level
#     }
#     RVAL <- list(statistic = STATISTIC,
#                  parameter = PARAMETER,
#                  p.value = as.numeric(PVAL),
#                  estimate = ESTIMATE,
#                  null.value = NVAL,
#                  alternative = alternative,
#                  method = method,
#                  data.name = DNAME)
#     if(conf.int)
#       RVAL <- c(RVAL, list(conf.int = cint))
#     class(RVAL) <- "htest"
#     RVAL
#   } else {
#     stop("Only Support Pearson now")
#   }
#   
# }
