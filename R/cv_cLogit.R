######################################################################
## This function is a modified version of the cv.lognet function from
## the glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate
#   Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.

## This program is free software: you can redistribute it and/or 
## modify it under the terms of the GNU General Public License as 
## published by the Free Software Foundation, either version 2 of 
## the License, or (at your option) any later version.

## They are copied here and modified to handle regression with multiple 
## independent studies (multi-task) and group/(sparse group) penalisation 

cv_cLogit <-
  function (outlist, lambda, x, y, foldid, type.measure,keep = FALSE)
  {
    prob_min = 1e-05
    prob_max = 1 - prob_min
    ## -- Objective for output -- ##
    K <- length(x)
    N = sapply(x, function(i) dim(i)[1])
    
    y <- lapply(y, function(y) {
      y = as.factor(y)
      ntab = table(y)
      nc = as.integer(length(ntab))
      diag(nc)[as.numeric(y), ]})
    
    nfolds = max(sapply(foldid, max))
    if (any(N/nfolds < 10) && type.measure == "auc") {
      warning("Too few (< 10) observations per fold for type.measure='auc' in cv.lognet; changed to type.measure='deviance'. Alternatively, use smaller value for nfolds",
              call. = FALSE)
      type.measure = "deviance"
    }
    nlams = length(lambda)
    predmat = matrix(NA, sum(N), nlams)
    foldidc <- do.call("c", foldid)
    for (i in seq(nfolds)) {
      which = lapply(1:K, function(k) foldid[[k]] == i)
      x_sub = lapply(1:K, function(k) x[[k]][which[[k]], , drop = FALSE])
      
      fitobj = outlist[[i]]
      preds = predict(fitobj,x_sub, type = "response")
      whichi = foldidc == i
      predmat[whichi, 1:nlams] = do.call('cbind',preds)
    }
    if (type.measure == "auc") {
      cvraw = matrix(NA, nfolds, nlams)
      for (i in seq(nfolds)) {
        which = foldidc == i
        for (j in seq(nlams)) {
          cvraw[i, j] = glmnet::auc.mat(y[which, ], predmat[which, j])
        }
      }
      N = nfolds
    }
    else {
      y <- do.call("rbind", y)
      ywt = apply(y, 1, sum)
      y = y/ywt
      N = nrow(y) - apply(is.na(predmat), 2, sum)
      cvraw = switch(type.measure,
                     mse = (y[, 1] - (1 - predmat))^2 +
                       (y[, 2] - predmat)^2,
                     mae = abs(y[, 1] - (1 - predmat)) +
                       abs(y[, 2] - predmat),
                     deviance = {
                       predmat = pmin(pmax(predmat, prob_min), prob_max)
                       lp = y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
                       ly = log(y)
                       ly[y == 0] = 0
                       ly = drop((y * ly) %*% c(1, 1))
                       2 * (ly - lp)
                     },
                     class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <=0.5)
      )
    }
    cvm = apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean,
                      na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, type.measure=type.measure)
    if (keep)
      out$fit.preval = predmat
    out
  }
