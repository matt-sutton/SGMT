######################################################################
## These functions are modified versions of the cv.glmnet function from
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

## The functions are copied here from glmnet and modified to 
## handle tuning over the alpha parameter in the group/sparse group selection. 


cv_alphas <- function(x, y, alphas, lambda = 0.01, verbose = F, nfolds = 10, nlam = 50, ...){
  K = length(x);  N = lapply(x, function(i) dim(i)[1])
  foldid = lapply(1:K, function(k) sample(rep(seq(nfolds), length = N[[k]])))
  betas = as.list(seq(alphas))
  cvsd<- cvm <- cvup <- cvlo <- matrix(0, nrow = nlam, ncol = length(alphas))
  lambdas = matrix(0, nrow = nlam, ncol = length(alphas))

  for(i in seq(alphas)){
    if(verbose) cat("tuning alpha =", alphas[i], "\n")
    temp <- cv_multi(x, y, lambda = lambda, foldid = foldid,
                     alpha = alphas[i], nlam=nlam, ...)
    cvm[,i] <- temp$cvm
    cvlo[,i] <- temp$cvlo
    cvup[,i] <- temp$cvup
    cvsd[,i] <- temp$cvsd
    lambdas[,i] <- temp$lambda
    betas[[i]] <- temp$multi.fit$beta
  }
  lamin=get_min_al(alphas, lambdas, cvm, cvsd)
  alpha.min = which(alphas==lamin$alpha.min)
  lambda.min = which(lambdas[,alpha.min]==lamin$lambda.min)
  lambda.1se = which(lambdas[,alpha.min]==lamin$lambda.1se)

  beta_min = betas[[lamin$alpha_min_ind[1]]][,,lamin$lambda_min_ind[1]]
  beta_1se = betas[[lamin$alpha_min_ind[1]]][,,lamin$lambda_1se_ind[1]]

  res <- list(cvm = cvm, cvsd=cvsd, cvlo = cvlo, cvup=cvup, lambdas=lambdas,
              alphas=alphas, beta_min=beta_min,beta_1se=beta_1se, betas = betas)
  res <- c(res, as.list(lamin))
  return(res)
}
get_min_al <-
  function(alpha,lambda,cvm,cvsd){
    if(is.null(dim(cvm))){
      lambda <- matrix(lambda, ncol = 1)
      cvm <- matrix(cvm, ncol = 1)
      cvsd <- matrix(cvsd, ncol = 1)
    }
    cv_min=min(cvm,na.rm=TRUE)
    idmin=cvm<=cv_min
    alpha_min_ind <- which(cvm<=cv_min, arr.ind = T)[,2]
    lambda_min_ind <- lambda_1se_ind <- rep(0, length(alpha_min_ind))
    lambda_min <- lambda_1se <- rep(0, length(alpha_min_ind))
    for( i in 1:length(alpha_min_ind)){
      inds_min <- cvm[,alpha_min_ind[i]]<=cv_min
      lambda_min[i] <- max(lambda[inds_min, alpha_min_ind[i]])
      lambda_min_ind[i] <- which(lambda_min[i] == lambda[, alpha_min_ind[i]])
      se_min = cvm[inds_min,alpha_min_ind[i]]+cvsd[inds_min,alpha_min_ind[i]]
      
      inds_se <- cvm[,alpha_min_ind[i]] <= se_min
      lambda_1se[i] <- max(lambda[inds_se, alpha_min_ind[i]])
      lambda_1se_ind[i] <- which(lambda_1se[i] == lambda[, alpha_min_ind[i]])
    }
    if( length(unique(alpha)) == 1 ){
      return(list(lambda.min=lambda_min[1],lambda.1se=lambda_1se[1]))
    } else {
      return(list(alpha_min_ind = alpha_min_ind, lambda_min_ind = lambda_min_ind, lambda_1se_ind = lambda_1se_ind,
                  alpha.min=alpha[alpha_min_ind],lambda.min=lambda_min,lambda.1se=lambda_1se))
    }
  }
getmin <- 
  function(lambda,cvm,cvsd){
    return(get_min_al(1,lambda,cvm,cvsd))
  }
