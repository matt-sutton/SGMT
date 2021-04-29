cLogit <- function(X, Y, lambda=0.0001, alpha = 1, groups, grp_weights = sqrt(table(groups)),
                   ind_weights = rep(1,length(groups)), rho = 1, abstol = 1e-04, reltol = 1e-03,
                   maxI = 10^2, lambdarel = 0.9, intercept = T, nlam=50, verbose = F, alpha_ADMM = 1){

  get_lambda <- function(lambda.mult){
    Xs <- lapply(X, function(x) x/nrow(x))
    X <- Matrix::bdiag(Xs)
    Y <- matrix(do.call("c", Y), ncol = 1)
    groups <- rep(groups,K)
    # ind_weights <- rep(ind_weights,K)
    # grp_weights <- rep(grp_weights,K)
    inds <- rep(1:p, K)
    lmaxlas <- sapply(1:p, function(i) norm(t(Y - Y*(1-Y)) %*% X[, i == inds]/ind_weights[i], "2"))
    lmaxgroup <- sapply(groups, function(i) norm(t(Y - Y*(1-Y)) %*% X[, i == groups]/grp_weights[i], "2"))

    ## Find the lambda which sets all variables to zero
    lmax <- max(alpha*lmaxlas + (1-alpha)*lmaxgroup)
    ## The calculation can sometimes be too low so return a slightly reduced lambda max
    lmax <- lmax*lambda.mult
    return(lmax)
  }

  ## -- Objective for output -- ##
  K <- length(X)
  N = sapply(X, function(i) dim(i)[1]); q <- K
  p = dim(X[[1]])[2];

  if(length(lambda) == 1) {
    lambda_path <- log(seq(exp(get_lambda(lambdarel)),exp(lambda), length.out = nlam))
  }  else {
    lambda_path <- lambda
    nlam <- length(lambda)
  }

  Y <- lapply(Y, function(a) as.matrix(model.matrix(~0+y, data = data.frame(y = factor(a))))[,-1, drop = F])
  if(intercept){
    X <- lapply(X, function(x) cbind(1, x))
  }

  b <- matrix(0, nrow = p + intercept, ncol = q)
  Beta <- array(0, dim = c(p+intercept, q, nlam))
  sol_info <- data.frame(rnorm = rep(0,nlam), snorm = rep(0,nlam), iter = rep(0,nlam))

  for (lam in 1:nlam){
    lambda <- lambda_path[lam]
    if(verbose){
      cat("lambda num:",lam,"\n")
    }
    z <- u <- matrix(0, nrow = p + intercept, ncol = q)
    for (i in 1:maxI){

      b <- lapply(1:K, function(k){
        bup <- b_update(X[[k]], Y[[k]], u[,k,drop=F] - z[,k, drop=F], rho, b[,k,drop = F])
        if(bup$converged) bup$coefficients else 0*b[,k,drop=F]
      })
      b <- do.call(cbind, b)

      # Update z with relaxation
      zold <- z
      bhat <- alpha_ADMM*b + (1-alpha_ADMM)*zold
      z <- weighted_proxop(x = bhat + u, groups = groups, lambda = mean(N)*lambda/rho, alpha = alpha,
                           grp_weights=grp_weights, ind_weights=ind_weights)

      u <- u + bhat - z

      #-- terminal checks (taken from Boyd matlab example)--#
      history_r_norm <- norm(b - z,"f");
      history_eps_pri <- sqrt(p)*abstol + reltol*max(norm(b), norm(-z))
      history_s_norm <- norm(rho*(z - zold));
      history_eps_dual <- sqrt(p)*abstol + reltol*norm(u)

      if(history_r_norm < history_eps_pri &&
         history_s_norm < history_eps_dual ){
        break
      }

      # Update penalty (Boyd 2010)
      if(history_r_norm > 10*history_s_norm){
        rho <- rho * 2
        u <- u/2
      }
      if(history_s_norm > 10*history_r_norm){
        rho <- rho / 2
        u <- u*2
      }
    }
    Beta[ , , lam] <- z
    sol_info[lam,] <-  c(history_r_norm, history_s_norm, i)
  }

  df <- apply(Beta, 3, function(a) sum(rowSums(abs(a)>0))) - intercept

  res <- list(beta = Beta,
              lambdas = lambda_path,
              sol_info = sol_info,
              df = df,
              intercept = intercept)
  class(res) = "cLogit"
  return(res)
}
