######################################################################
## This function is a modified version of the cv.glmnet function from
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
## This function calls cv_cLogit when performing cross validation for a set 
## alpha lambda combination.

cv_multi <-
  function (x, y, lambda = 0.0001, groups = rep(1,ncol(x[[1]])), type.measure = c("mse","deviance", "class", "auc", "mae"),
            nfolds = 10, foldid, keep = FALSE, parallel = FALSE, verbose = F, maxI = 500, ...)
  {
    if (missing(type.measure))
      type.measure = "deviance"
    else type.measure = match.arg(type.measure)
    #-- params --#
    K = length(x)
    N = lapply(x, function(i) dim(i)[1]);
    #-- Call --#
    multi.call = match.call(expand.dots = TRUE)
    which = match(c("type.measure", "nfolds", "foldid", "grouped",
                    "keep"), names(multi.call), F)
    if (any(which))
      multi.call = multi.call[-which]
    multi.call[[1]] = as.name("cLogit")
    start_time <- Sys.time()
    multi.obj = cLogit(x, y, lambda = lambda, groups = groups, maxI = maxI, ...)
    end_time <- Sys.time()
    multi.obj$call = multi.call
    timefit1 <- end_time - start_time

    subclass=class(multi.obj)[[1]]
    lambda=multi.obj$lambda
    df = multi.obj$df
    if (missing(foldid))
      foldid = lapply(1:K, function(k) sample(rep(seq(nfolds), length = N[[k]])))
    else nfolds = max(sapply(foldid, max))
    if(verbose){
      cat("\nTime for single fit ", timefit1, "Est time: ", timefit1*nfolds," \n")
    }
    if (nfolds < 3)
      stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist = as.list(seq(nfolds))
    if (parallel) {
      #  if (parallel && require(foreach)) {
      outlist = foreach(i = seq(nfolds), .packages = c("pleioMethods")) %dopar%
        {
          which = lapply(1:K, function(k) foldid[[k]] == i)
          y_sub = lapply(1:K, function(k) y[[k]][!which[[k]], , drop = FALSE])
          x_sub = lapply(1:K, function(k) x[[k]][!which[[k]], , drop = FALSE])
          cLogit(x_sub, y_sub, lambda = lambda, groups = groups, maxI = maxI, ...)
        }
    }
    else {
      for (i in seq(nfolds)) {
        if(verbose){
          cat(" fold",i)
        }
        which = lapply(1:K, function(k) foldid[[k]] == i)
        y_sub = lapply(1:K, function(k) y[[k]][!which[[k]], , drop = FALSE])
        x_sub = lapply(1:K, function(k) x[[k]][!which[[k]], , drop = FALSE])
        outlist[[i]] = cLogit(x_sub, y_sub, lambda = lambda, groups = groups, maxI = maxI, ...)
      }
    }
    fun = paste("cv", subclass, sep = "_")
    lambda = multi.obj$lambda
    cvstuff = do.call(fun, list(outlist, lambda, x, y, foldid, type.measure, keep))
    cvm = cvstuff$cvm
    cvsd = cvstuff$cvsd
    nas=is.na(cvsd)
    if(any(nas)){
      lambda=lambda[!nas]
      cvm=cvm[!nas]
      cvsd=cvsd[!nas]
      df=df[!nas]
    }
    cvname = cvstuff$type.measure
    out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
                 cvsd, cvlo = cvm - cvsd, name = cvname, multi.fit = multi.obj, nzero=df)
    if (keep)
      out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
    lamin=if(cvname=="AUC")getmin(lambda,-cvm,cvsd)
    else getmin(lambda, cvm, cvsd)
    whichlambdamin = lapply(as.list(lamin), function(l) which(l == lambda))
    names(whichlambdamin) = c("lambda.mini", "lambda.1sei")
    obj = c(out, as.list(lamin), whichlambdamin)
    class(obj) = "cv.multi"
    obj
  }
