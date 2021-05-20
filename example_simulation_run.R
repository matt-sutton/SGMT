############################################################
#                                                          #
####       Code replication example                     ####
# This code replicates the analysis for a run of the       #
# simulation study. Here we show a plots of the regression #
# recovery performance for our method compared with ASSET  #
# The script demonstrates how bootstrapping was used to    #
# select variables for the methods and how FDR was set for #
# ASSET                                                    #
#                                                          #
##%######################################################%##


library(parallel)
library(foreach)
library(tidyverse)
library(SGMT)

##%######################################################%##
#                                                          #
####                   Help function                    ####
#                                                          #
##%######################################################%##

bootstrap_data = function(X, Y) {
  for( i in 1:length(X)){
    nobs = nrow(X[[i]])
    newIndex = sample(seq(nobs), size = nobs, replace = TRUE)
    X[[i]] = X[[i]][newIndex, ]
    Y[[i]] = Y[[i]][newIndex, ]
  }
  return(list(X, Y))
}
make_summary_stat = function(X, Y, covar=NULL){
  require(tidyverse)
  
  # Count non missing case and control for each SNP
  case = apply(X[Y==1,], 2, function(x) length(!is.na(x)))
  control = apply(X[Y==0,], 2, function(x) length(!is.na(x)))
  n_case_control = enframe(case, name="snp", value="ncase") %>%
    left_join(enframe(control, name="snp", value="ncontrol"), by=c("snp")) %>%
    mutate(n = ncase + ncontrol)
  
  # Run logistic regression on each SNP with apply
  t0 = Sys.time()
  
  tmp = apply(X, 2, function(x) {
    x = matrix(x, ncol = 1)
    if (is.null(covar)) {
      fit = glm(Y ~ x, family=binomial())
    } else {
      fit = glm(Y ~ covar + x, family=binomial())
    }
    co = summary(fit)$coefficients
    return(c(co["x", "Estimate"], co["x", "Std. Error"], co["x", "Pr(>|z|)"]))
  })
  
  summary_stat = t(tmp) %>% as_tibble(rownames = "snp") %>%
    rename(beta=V1, sigma=V2, pval=V3) %>%
    mutate(z = beta / sigma) %>%
    left_join(n_case_control, by=c("snp"))
  
  Sys.time() - t0
  
  return(summary_stat)
}

make_summary_stat2 = function(X, Y, covar=NULL){
  require(tidyverse)
  
  # Run logistic regression on each SNP with apply
  t0 = Sys.time()
  
  tmp = apply(X, 2, function(x) {
    x = matrix(x, ncol = 1)
    if (is.null(covar)) {
      fit = glm(Y ~ x, family=binomial())
    } else {
      fit = glm(Y ~ covar + x, family=binomial())
    }
    co = summary(fit)$coefficients
    return(c(co["x", "Estimate"], co["x", "Std. Error"], co["x", "Pr(>|z|)"]))
  })
  
  summary_stat = t(tmp) %>% as_tibble(rownames = "snp") %>%
    rename(beta=V1, sigma=V2, pval=V3) %>%
    mutate(z = beta / sigma)
  
  
  Sys.time() - t0
  
  return(summary_stat)
}

##%######################################################%##
#                                                          #
####                1.2 Plot functions                ####
#                                                          #
##%######################################################%##

plot_beta = function(beta, groups, main = "Beta") {
  for(s in 1:ncol(beta)){
    plot(beta[,s], col=groups, type="h", main=paste0(main," - Study ",s), ylab = '')
  }
}

##%######################################################%##
#                                                          #
####                     Simulation                     ####
#                                                          #
##%######################################################%##


##%######################################################%##
# Simulation setting

library(doParallel)
registerDoParallel(cores = 6)
getDoParWorkers()

data <- gen_data_paper(scenario = 2, seed = 1)
grpvec <- data$groups
grp_size <- table(grpvec)
true_beta <- data$beta

par(mfrow = c(1,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
par(mfrow = c(1,1))

X <- data$X
Y <- data$Y
P <- ncol(X[[1]])
K <- length(X)
ngrp <- length(unique(grpvec))
methods = list()

## You can use different choices for the weights. Eg the ols solution
var_weights <- c()
for(i in 1:P){
  study_w <- c()
  for(j in 1:K){
    ## get ols from each var ~~~
    Xtemp <- X[[j]][,i]
    beta_ols <- coef(glm(y~0+., family = "binomial", data = data.frame(y = Y[[j]],x = Xtemp)))
    study_w <- c(study_w,beta_ols)
  }
  var_weights <- c(var_weights,1/norm(matrix(study_w),type="2"))
}
# Different group weights are possible too.
grp_weights=rep(NA,ngrp)
for(j in 1:ngrp){
  grp_weights[j]=1/sqrt(sum((1/var_weights[which(grpvec==j)])**2))
}

methods$data$trueBeta = true_beta
methods$data$X = X
methods$data$Y = Y
methods$data$grpvec = grpvec

cat("\nStart SMT \n")
## Boot approach
t0 = Sys.time()
resSMT = cv_multi(X, Y, alpha = 1, groups = grpvec, intercept = F, verbose=T, ind_weights=var_weights, grp_weights = grp_weights)
t1 = Sys.time()

delta_time = difftime(t1, t0, units = "mins")
methods$SMT$time = delta_time
methods$SMT$lambda_min = resSMT$lambda.min
methods$SMT$beta_min = resSMT$multi.fit$beta[,,resSMT$lambda.mini]
methods$SMT$lambda_1se = resSMT$lambda.1se
methods$SMT$beta_1se = resSMT$multi.fit$beta[,,resSMT$lambda.1sei]

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(methods$SMT$beta_1se, groups = grpvec, main = "SMT")
par(mfrow = c(2,1))

#######################################
## Boot ##
bootR <- foreach ( b = 1:200, .packages = c("parallel","tidyverse","SGMT"), .combine = cbind ) %dopar% {
  b_XY <- bootstrap_data(X, Y)
  ## SMT
  resBootStrap = cLogit(b_XY[[1]], b_XY[[2]], ind_weights = var_weights, grp_weights = grp_weights,
                        lambda=c(methods$SMT$lambda_1se, 0), alpha = 1,
                        groups = grpvec, intercept = FALSE)
  # Extract beta and save
  abs(c(resBootStrap$beta[,,1])) > 1e-12 ## Threshold from
}
methods$SMT$Inc_prob <- rowMeans(bootR)

## Plot selected using bootstrap (strategy 1)
inc_var <- which(matrix(methods$SMT$Inc_prob == 1, P, 2), arr.ind = T)
SMT_beta_est_1 <- methods$SMT$beta_1se*0
SMT_beta_est_1[inc_var] <- methods$SMT$beta_1se[inc_var]

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(SMT_beta_est_1, groups = grpvec, main = "SMT (BS1)")
par(mfrow = c(1,1))

## Plot selected using bootstrap (strategy 2)
non_selected <- which(abs(methods$SMT$beta_1se) < 1e-13)
acc_rate <- max(methods$SMT$Inc_prob[non_selected]) # Max selection rate of non-selected
inc_var_app <- which(matrix(methods$SMT$Inc_prob > acc_rate, P, 2), arr.ind = T)
SMT_beta_est_2 <- methods$SMT$beta_1se*0
SMT_beta_est_2[inc_var_app] <- methods$SMT$beta_1se[inc_var_app]

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(SMT_beta_est_2, groups = grpvec, main = "SMT (BS2)")
par(mfrow = c(1,1))

#### 2.2 GMT ####
cat("\nStart GMT \n")
t0 = Sys.time()
resGMT = cv_multi(X, Y, alpha = 0, groups = grpvec, intercept = F, verbose=T, ind_weights=var_weights, grp_weights = grp_weights)
t1 = Sys.time()

delta_time = difftime(t1, t0, units = "mins")
methods$GMT$time = delta_time
methods$GMT$lambda_min = resGMT$lambda.min
methods$GMT$beta_min = resGMT$multi.fit$beta[,,resGMT$lambda.mini]
methods$GMT$lambda_1se = resGMT$lambda.1se
methods$GMT$beta_1se = resGMT$multi.fit$beta[,,resGMT$lambda.1sei]

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(methods$GMT$beta_1se, groups = grpvec, main = "GMT")
par(mfrow = c(2,1))

bootR <- foreach ( b = 1:200, .packages = c("parallel","tidyverse","SGMT"), .combine = cbind ) %dopar% {
  b_XY <- bootstrap_data(X, Y)
  ## SMT
  resBootStrap = cLogit(b_XY[[1]], b_XY[[2]], ind_weights = var_weights, grp_weights = grp_weights,
                        lambda=c(methods$GMT$lambda_1se, 0), alpha = 0,
                        groups = grpvec, intercept = FALSE)
  # Extract beta and save
  abs(c(resBootStrap$beta[,,1])) > 1e-12
}
methods$GMT$Inc_prob <- rowMeans(bootR)

## Plot selected using bootstrap (strategy 1)
inc_var <- which(matrix(methods$GMT$Inc_prob == 1, P, 2), arr.ind = T)
GMT_beta_est_1 <- methods$GMT$beta_1se*0
GMT_beta_est_1[inc_var] <- methods$GMT$beta_1se[inc_var]

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(GMT_beta_est_1, groups = grpvec, main = "GMT (BS1)")
par(mfrow = c(1,1))

## Plot selected using bootstrap (strategy 2)
non_selected <- which(abs(methods$GMT$beta_1se) < 1e-13)
acc_rate <- max(methods$GMT$Inc_prob[non_selected]) # Max selection rate of non-selected
inc_var_app <- which(matrix(methods$GMT$Inc_prob > acc_rate, P, 2), arr.ind = T)
GMT_beta_est_2 <- methods$GMT$beta_1se*0
GMT_beta_est_2[inc_var_app] <- methods$GMT$beta_1se[inc_var_app]

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(GMT_beta_est_2, groups = grpvec, main = "GMT (BS2)")
par(mfrow = c(1,1))


#### 2.3 SGMT ####
cat("\nStart SGMT \n")
t0 = Sys.time()
resSGMT = cv_alphas(X, Y, alphas = seq(0.05,0.95, by = 0.1), groups = grpvec, intercept = F, verbose=T, ind_weights=var_weights, grp_weights = grp_weights)
t1 = Sys.time()
# compare_beta(resSGMT$beta_1se, true_beta, groups = grpvec)
delta_time = difftime(t1, t0, units = "mins")
methods$SGMT$time = delta_time
methods$SGMT$alpha_min = resSGMT$alpha.min
methods$SGMT$lambda_min = resSGMT$lambda.min
methods$SGMT$beta_min = resSGMT$beta_min
methods$SGMT$alpha_1se = resSGMT$alpha.min
methods$SGMT$lambda_1se = resSGMT$lambda.1se
methods$SGMT$beta_1se = resSGMT$beta_1se

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(methods$SGMT$beta_1se, groups = grpvec, main = "SGMT")
par(mfrow = c(2,1))

bootR <- foreach ( b = 1:200, .packages = c("parallel","tidyverse","SGMT"), .combine = cbind ) %dopar% {
  b_XY <- bootstrap_data(X, Y)
  ## SMT
  resBootStrap = cLogit(b_XY[[1]], b_XY[[2]], ind_weights = var_weights, grp_weights = grp_weights,
                        lambda=c(methods$SGMT$lambda_1se, 0), alpha = methods$SGMT$alpha_1se,
                        groups = grpvec, intercept = FALSE)
  # Extract beta and save
  abs(c(resBootStrap$beta[,,1])) > 1e-12
}
methods$SGMT$Inc_prob <- rowMeans(bootR)

## Plot selected using bootstrap (strategy 1)
inc_var <- which(matrix(methods$SGMT$Inc_prob == 1, P, 2), arr.ind = T)
SGMT_beta_est_1 <- methods$SGMT$beta_1se*0
SGMT_beta_est_1[inc_var] <- methods$SGMT$beta_1se[inc_var]

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(SGMT_beta_est_1, groups = grpvec, main = "SGMT (BS1)")
par(mfrow = c(1,1))

## Plot selected using bootstrap (strategy 2)
non_selected <- which(abs(methods$SGMT$beta_1se) < 1e-13)
acc_rate <- max(methods$SGMT$Inc_prob[non_selected]) # Max selection rate of non-selected
inc_var_app <- which(matrix(methods$SGMT$Inc_prob > acc_rate, P, 2), arr.ind = T)
SGMT_beta_est_2 <- methods$SGMT$beta_1se*0
SGMT_beta_est_2[inc_var_app] <- methods$SGMT$beta_1se[inc_var_app]

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(SGMT_beta_est_2, groups = grpvec, main = "SGMT (BS2)")
par(mfrow = c(1,1))




#### 2.7 ASSET ####
#BiocManager::install("ASSET")
cat("\nStart ASSET \n")
t0 = Sys.time()
summary_stat_S1 = make_summary_stat(X[[1]], Y[[1]])
summary_stat_S2 = make_summary_stat(X[[2]], Y[[2]])
beta = as.matrix(data.frame(Study1 = summary_stat_S1$beta,
                            Study2 = summary_stat_S2$beta,
                            row.names = summary_stat_S1$snp))
sigma = as.matrix(data.frame(Study1 = summary_stat_S1$sigma,
                             Study2 = summary_stat_S2$sigma,
                             row.names = summary_stat_S1$snp))
case = as.matrix(data.frame(Study1 = summary_stat_S1$ncase,
                            Study2 = summary_stat_S2$ncase,
                            row.names = summary_stat_S1$snp))
control = as.matrix(data.frame(Study1 = summary_stat_S1$ncontrol,
                               Study2 = summary_stat_S2$ncontrol,
                               row.names = summary_stat_S1$snp))
SNPs = summary_stat_S1$snp
Study = c("Study1", "Study2")

library(ASSET)
res_ASSET = h.traits(SNPs, Study, beta, sigma, case, control, meta=TRUE)
res_2sides = h.summary(res_ASSET)$Subset.2sided

t1 = Sys.time()
delta_time = difftime(t1, t0, units = "mins")
methods$ASSET$time = delta_time
methods$ASSET$raw_res = res_2sides

ASSET_res = methods$ASSET$raw_res

pleiotropy = rep(0,dim(ASSET_res)[1])
sig_p1=which(ASSET_res$Pvalue.1<0.05)
sig_p2=which(ASSET_res$Pvalue.2<0.05)

pleiotropy[sig_p1]=pleiotropy[sig_p1]+
  as.numeric(lapply(strsplit(as.character(ASSET_res[sig_p1,"Pheno.1"]),split=',',
                             fixed=TRUE),length))
pleiotropy[sig_p2]=pleiotropy[sig_p2]+
  as.numeric(lapply(strsplit(as.character(ASSET_res[sig_p2,"Pheno.2"]),split=',',
                             fixed=TRUE),length))
ASSET_edit=cbind(ASSET_res,pleiotropy)
colnames(ASSET_edit)<-c('snp','pval','pval_p','pval_n','OR_p','CI.l.p',
                        'CI.h.p','OR_n','CI.l.n','CI.h.n',
                        'pheno.p','pheno.n','pleiotropy')
annot <- data.frame(snp = ASSET_edit$snp, gene = grpvec)
res_annot_ASSET=merge(ASSET_edit, annot, by="snp", sort = FALSE)

## Not Accounting for multiple comparisons
estimBeta = res_annot_ASSET %>%
  mutate(pheno.p = ifelse(pval < 0.05, pheno.p, "")) %>%
  mutate(pheno.n = ifelse(pval < 0.05, pheno.n, "")) %>%
  mutate(beta_study1 = 0) %>%
  mutate(beta_study1=ifelse(str_detect(pheno.p,"Study1"),log(OR_p),beta_study1)) %>%
  mutate(beta_study1=ifelse(str_detect(pheno.n,"Study1"),log(OR_n),beta_study1)) %>%
  mutate(beta_study2 = 0) %>%
  mutate(beta_study2=ifelse(str_detect(pheno.p,"Study2"),log(OR_p),beta_study2)) %>%
  mutate(beta_study2=ifelse(str_detect(pheno.n,"Study2"),log(OR_n),beta_study2)) %>%
  dplyr::select(beta_study1, beta_study2) %>% as.matrix() %>% c()

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(matrix(estimBeta,P,2), groups = grpvec, main = "ASSET")
par(mfrow = c(1,1))

#####
# FDR by SNPs
FDR_SNP=p.adjust(res_annot_ASSET$pval, method = 'fdr')
res_annot_ASSET=cbind(res_annot_ASSET,FDR_SNP)

estimBeta_FDR = res_annot_ASSET %>%
  mutate(pheno.p = ifelse(FDR_SNP < 0.05, pheno.p, "")) %>%
  mutate(pheno.n = ifelse(FDR_SNP < 0.05, pheno.n, "")) %>%
  mutate(beta_study1 = 0) %>%
  mutate(beta_study1=ifelse(str_detect(pheno.p,"Study1"),log(OR_p),beta_study1)) %>%
  mutate(beta_study1=ifelse(str_detect(pheno.n,"Study1"),log(OR_n),beta_study1)) %>%
  mutate(beta_study2 = 0) %>%
  mutate(beta_study2=ifelse(str_detect(pheno.p,"Study2"),log(OR_p),beta_study2)) %>%
  mutate(beta_study2=ifelse(str_detect(pheno.n,"Study2"),log(OR_n),beta_study2)) %>%
  dplyr::select(beta_study1, beta_study2) %>% as.matrix() %>% c()

par(mfrow = c(2,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
plot_beta(matrix(estimBeta_FDR,P,2), groups = grpvec, main = "ASSET (Selected)")
par(mfrow = c(1,1))

## Compare All Methods --- Bootstrap sampling strategy 1 (simulation)
par(mfrow = c(5,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
## SGMT
plot_beta(SGMT_beta_est_1, groups = grpvec, main = "SGMT (BS1)")
## GMT
plot_beta(GMT_beta_est_1, groups = grpvec, main = "GMT (BS1)")
## SMT
plot_beta(SMT_beta_est_1, groups = grpvec, main = "SMT (BS1)")
## ASSET
plot_beta(matrix(estimBeta_FDR,P,2), groups = grpvec, main = "ASSET (FDR)")
par(mfrow = c(1,1))

## Compare All Methods  --- Bootstrap sampling strategy 2 (application)
par(mfrow = c(5,2), mar = c(3,2,2,2))
plot_beta(true_beta, groups = grpvec)
## SGMT
plot_beta(SGMT_beta_est_2, groups = grpvec, main = "SGMT (BS2)")
## GMT
plot_beta(GMT_beta_est_2, groups = grpvec, main = "GMT (BS2)")
## SMT
plot_beta(SMT_beta_est_2, groups = grpvec, main = "SMT (BS2)")
## ASSET
plot_beta(matrix(estimBeta_FDR,P,2), groups = grpvec, main = "ASSET (FDR)")
par(mfrow = c(1,1))


#save.image("example_run.Rdata")

library(randomcoloR)
n <- length(resSGMT$alphas)
palette <- distinctColorPalette(n)
par(mfrow = c(1,1), mar = rep(4,4))
m <- matrix(c(1,2,1,3),2,2)
layout(m)

plot_cv(resSGMT, cols = palette)
title("SGMT")
plot_cv(resSMT)
title("SMT")
plot_cv(resGMT)
title("GMT")