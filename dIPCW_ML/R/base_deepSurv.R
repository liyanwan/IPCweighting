#' A helper function
getCost <- function (y, yh, loss = "mse", weights = NULL) {
  n = length(yh)
  if (loss == "mse") {
    if (!is.null(weights)) {
      if (length(weights) != n) 
        stop("Error: Wrong size for weights")
    }
    else weights = 1
    dy = weights * (y - yh)
    cost = sum(weights * (y - yh)^2)/2
  }
  else if (loss == "cox") {
    delta = y[, 2]
    eb = exp(yh)
    S0 = cumsum(eb)
    ht = delta/S0
    Ht = .rcumsum(ht)
    cost = -sum(delta * (yh - log(S0)))
    dy = (delta - eb * Ht)
  }
  else if (loss == "log") {
    y = y * weights
    yh = yh * weights
    cost = sum(yh - y * log(yh))
    dy = (yh - y)/yh
  }
  else if (loss == "bin") {
    eps <- 1e-5
    yh   <- pmin(pmax(yh, eps), 1 - eps)  
    cost = sum((y * log(yh) + (1 - y) * log(1 - yh)) * weights)
    dy = (y/yh - (1 - y)/(1 - yh)) * weights
  }
  else if (loss == "mae") {
    dy = weights * abs(y - yh)
    cost = sum(dy)
    dy = ifelse(dy > 0, -1, 1)
  }
  else {
    stop(paste("loss function", loss, "is not defiend yet"))
  }
  return(list(cost = cost, dy = dy))
}

#' K-fold CV for deepSurv using C-index or NPLL (partial log-likelihood)
#'
#' @param X numeric data.frame or matrix of predictors (n x p)
#' @param Y data.frame with columns: `observed_time` (time), `delta` (status 1/0), and optional `E` (ignored).
#' @param alpha numeric vector (momentum/optimizer param)
#' @param lambda numeric vector (L2 weight decay)
#' @param lr_rate numeric vector (learning rate)
#' @param model list; dnn model definition as expected by `dnn::dNNmodel`
#' @param measure one of c("cindex", "npll")
#' @param k The number of folds used in cross-validation.
#' @param foldid A vector specifying the fold assignment for each observation.
#' @param verbose logical; print training progress/cost messages from dnn
#' @param sample_batch If TRUE, ensures that the same set of samples is selected in each batch across cross-validation folds.
#' @param batch_seed Used with set.seed(batch_seed) to guarantee reproducibility of batch sampling when sample_batch = TRUE.

CV_deepSurv <- function(X, 
                        Y, 
                        alpha,
                        lambda,
                        lr_rate,
                        model,
                        measure,
                        k = 5, 
                        foldid = NULL,
                        epochs = 200,
                        batch_size = 64,
                        epsilon = 1e-3,
                        verbose = FALSE,
                        sample_batch = TRUE,
                        batch_seed = 0,
                        ...) { 
  if (is.data.frame(X)) {
    if (!all(vapply(X, is.numeric, logical(1)))) stop("All columns of X must be numeric.")
    X = as.matrix(X)
  } else if (!is.matrix(X)) {
    X = as.matrix(X)
  }
  if (is.data.frame(Y)) {
    req = c("observed_time", "sigma")
    if (!all(req %in% colnames(Y))) stop("Y must have columns: observed_time, sigma")
    Y = as.matrix(Y[, req, drop = FALSE])
  } else if (is.matrix(Y)) {
    if (all(c("observed_time","delta") %in% colnames(Y))) {
      Y = as.matrix(Y[, c("observed_time","delta"), drop = FALSE])
    } else {
      Y = as.matrix(Y[, 1:2, drop = FALSE])
    }
  } else {
    stop("Y must be a data.frame or matrix")
  }
  stopifnot(nrow(X) == nrow(Y), ncol(Y) == 2)
  data = cbind(X, Y)
  var_name = colnames(X)
  cv_results = expand.grid(
    alpha = alpha,
    lambda = lambda,
    lr_rate = lr_rate
  )
  cv_results[[measure]] = NA
  
  # Folds
  if (is.null(foldid)) {
    foldid = sample(rep(1:k, length.out = nrow(X)))
  }
  k = max(foldid)
  fold_sizes = as.integer(table(factor(foldid, levels = seq_len(k))))
  fold_weights = fold_sizes / sum(fold_sizes)
  
  is_cindex = m %in% c("cindex")
  is_deviance = !is_cindex  # treat everything else as deviance/NPLL
  getCost = try(getFromNamespace(".getCost", "dnn"), silent = TRUE)
  have_internal_cost = !inherits(getCost, "try-error")
  for (i in seq_len(nrow(cv_results))) {
    print(i)
    al = cv_results$alpha[i]
    lam = cv_results$lambda[i]
    lr = cv_results$lr_rate[i]

    fold_performance_measure = rep(NA, k)
    for (f in seq_len(k)) {
      idx_tr = which(foldid != f)
      idx_te = which(foldid == f)
      train_data = data[idx_tr, , drop = FALSE]
      test_data  = data[idx_te, , drop = FALSE]
      train_X = as.matrix(train_data[, colnames(X), drop = FALSE])
      test_X  = as.matrix(test_data[, colnames(X), drop = FALSE])
      # Y assumed as last 2 columns of data (time,status); safer to slice by position
      y_cols = (ncol(data) - 1):ncol(data)
      train_Y = as.matrix(train_data[, y_cols, drop = FALSE])
      test_Y  = as.matrix(test_data[, y_cols, drop = FALSE])
      formula = as.formula(paste("Surv(observed_time, sigma) ~", paste(var_name, collapse = " + ")))
      if(sample_batch) set.seed(batch_seed)
      fit = dnn::deepSurv(formula, data = cbind.data.frame(train_X, train_Y), model = model, 
                           epochs = epochs, lr_rate = lr, batch_size = batch_size, 
                           alpha = al, lambda = lam, verbose = verbose)
      if (inherits(fit, "try-error")) next
      if (is_cindex) {
        val = predict(fit, test_X, test_Y)$c.index
      } else {
        lp = tryCatch({as.numeric(predict(fit, newdata = test_X)$predictors)},error  = function(e) {NA})
        if(all(is.na(lp))) {val = NA} else {
          ord = order(test_Y[, 1], decreasing = TRUE)
          y2 = test_Y[ord, , drop = FALSE]
          lp2 = lp[ord]
          val = getCost(y2, lp2, loss = "cox", weights = NULL)$cost
          }
        }
      if (!inherits(val, "try-error") && is.finite(val)) fold_performance_measure[f] = val
    }
    ok = is.finite(fold_performance_measure)
    if (all(ok)) {
      cv_results[[measure]][i] = sum(fold_performance_measure[ok] * w)
    }
  }
  if (is_cindex) {
    ord = order(-cv_results[[measure]])
  } else {
    ord = order(cv_results[[measure]])
  }
  cv_results = cv_results[ord, , drop = FALSE]
  return(list(results = cv_results, measure = measure, control_default = ctrl_defaults))
}



#' Grid search tuner for deepSurv based on CV (C-index or NPLL)

tune_deepSurv <- function(
    X, Y, X_validation = NULL,
    alpha, lambda, lr_rate,
    model,
    time_point,
    measure = c("cindex", "npll"),
    k = 5, foldid = NULL,
    epochs = 200, batch_size = 64, epsilon = 1e-3, verbose = FALSE, sample_batch = TRUE, batch_seed = 0,...) {
  measure = match.arg(measure)
  cv = CV_deepSurv(X, Y, alpha, lambda, lr_rate, model, measure, k, foldid,
                    epochs = epochs, batch_size = batch_size, epsilon = epsilon, sample_batch = sample_batch,
                    batch_seed = batch_seed, verbose = verbose)
  if(length(unique(cv$results[[measure]])) == 1){
    best = as.data.frame(t(vapply(cv$results,
                                   function(x) median(x, na.rm = TRUE),
                                   numeric(1))))
    } else {
    best = cv$results[1L, , drop = FALSE]
    }
  fit = NULL
  var_name = colnames(X)
  if (inherits(Y, "Surv")) Y = cbind(time = Y[, 1], status = Y[, 2])
  if (is.data.frame(Y)) Y = as.matrix(Y[,c(1,2)])
  formula = as.formula(paste("Surv(observed_time, sigma) ~", paste(var_name, collapse = " + ")))
  if(sample_batch) set.seed(batch_seed)
  fit = dnn::deepSurv(formula, data = cbind.data.frame(Y, X), model = model,
                       lr_rate = best$lr_rate, alpha = best$alpha, lambda = best$lambda,
                       epochs = epochs, batch_size = batch_size, epsilon = epsilon)
  if(!is.null(X_validation)){
    surv = survfit(Surv(Y[,'observed_time'], Y[,'sigma'])~1)
    t0 = summary(surv, times = time_point)
    base = t0$surv
    risk = predict(fit, newdata = as.matrix(X_validation))$risk
    coxph_penal_event_probability = 1-base^risk
  } else {
    coxph_penal_event_probability = NULL
  }
  return(data.frame(DeepSurv = coxph_penal_event_probability))
}



