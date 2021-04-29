invlogit <- function (x) {
  1/(1 + exp(-x))
}
sim_contX = function(n, p, sigma) {
  require(MASS)
  X = mvrnorm(n=n, mu=rep(0, p), Sigma=sigma)
  colnames(X) = paste0("rs", 1:p)
  return(X)
}
grp_list2vec = function(grp) {
  p = length(unlist(grp))
  grpv = vector(mode="integer", length = p)
  i=1
  for (g in grp) {
    grpv[g] = i
    i = i+1
  }
  return(grpv)
}
## Function to simulate data
gen_data_MultiLogit <- function( nsample = rep(60,5), beta = NULL){
  K <- length(nsample)
  X = lapply(1:K,function(k) matrix(rnorm(nsample[k]*p),nsample[k],p))
  Y = lapply(1:K, function(i) matrix(rbinom(nsample[i],1,invlogit(X[[i]]%*% beta[,i]))))
  list(X=X, Y=Y, beta=beta)
}

## Function to simulate data from the paper
gen_data_paper <- function( scenario = 2, seed = 1 ){
  
  scenario <- scenario*2 
  K <- 2;  N <- rep(scenario*50, K)
  P <- scenario*80;  ngrp <- scenario*4
  corr_in_grp <- 0; MAF <- 0.25
  n_true_var <- 2*scenario;  beta_order <- c(0.8,0.8)
  grp <- parallel::splitIndices(P, ngrp)
  grpvec <- grp_list2vec(grp)
  grp_size <- table(grpvec)
  
  # Make beta
  tbeta <- matrix(0,nrow=P,ncol=K)
  startInd <- 1
  for( i in 1:scenario){
    endInd <- startInd + 20 - 1
    tbeta[startInd:endInd,1] <- c(rep(beta_order[1],n_true_var),
                                  rep(0,grp_size[1]-n_true_var))
    tbeta[startInd:endInd,2] <- (-1)^(i+1)*c(rep(beta_order[2],n_true_var),
                                             rep(0,grp_size[1]-n_true_var))
    startInd <- endInd + 1
  }
  set.seed(seed)
  # Make Data
  X <- lapply(1:K, function(x) {
    sim_contX(N[x], P, diag(1,P,P))
  })
  Y <- lapply(1:K, function(x) {
    as.matrix(rbinom(N[x], 1, invlogit(X[[x]] %*% tbeta[,x])),ncol = 1)
  })
  data <- list(X=X, Y=Y, beta = tbeta, groups = grpvec)
  return(data)
}

