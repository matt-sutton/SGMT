predict.cLogit <- function(object, newx, lambda = 1:length(object$lambdas), type = "response"){
  
  nlam <- length(lambda)
  K <- length(newx)
  if(object$intercept){
    newx <- sapply(newx, function(a) cbind(1, a) )
  }
  fitProbs <- lapply(1:nlam, function(a) {
    xb <- lapply(1:K, function(k) newx[[k]]%*%object$beta[,k,a])
    xb <- do.call('rbind', xb)
    exp(xb)/(1+rowSums(exp(xb)))
  })
  fitVals <- lapply(fitProbs, function(a) (a < 0.5) + 0.0)
  if(type == "response"){
    return(fitProbs)
  } else {
    return(fitVals)
  }
}
