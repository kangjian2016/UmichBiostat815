# UmichBiostat815

##Overview

This package implements some statistical methods and algorithms that were taught in Biostatistics 815 at the University of Michgian in Winter 2023

##Install from Github

Install package `UmichBiostat815' from Github with 

```r
devtools::install_github("kangjian2016/UmichBiostat815")
library(UmichBiostat815)
```

##Examples: Ridge Regression

```r
library(glmnet)
set.seed(2023)
regdat <- simul_linear_reg(regcoef = c(rep(c(1,-1),length=10),rep(0,length=990)),
intercept = 2, 
n = 100, R_sq = 0.9)
lambda <- 0.04
res_glmnet <- with(regdat,glmnet(x=X,y=y,lambda=lambda,alpha=0.0))
res_solve <- with(regdat, linear_ridge_reg(y = y, X = X, lambda = lambda,method="solve"))
res_solve_p <- with(regdat, linear_ridge_reg(y = y, X = X, lambda = lambda,method="solve-large-p"))
loss <- list()
loss$glmnet <- with(res_glmnet,ridge_reg_loss(intercept=a0,regcoef=beta,lambda=lambda,regdat$y,regdat$X))
loss$solve <- with(res_solve,ridge_reg_loss(intercept,regcoef,lambda=lambda,regdat$y,regdat$X))
loss$solve_p <- with(res_solve_p,ridge_reg_loss(intercept,regcoef,lambda=lambda,regdat$y,regdat$X))
loss$true <- with(regdat,ridge_reg_loss(intercept,regcoef,lambda=lambda,y,X))
print(as.data.frame(loss))
mse <- list()
mse$glmnet <- mean((res_glmnet$beta - regdat$regcoef)^2)
mse$solve <- mean((res_solve$regcoef - regdat$regcoef)^2)
mse$solve_p <- mean((res_solve_p$regcoef - regdat$regcoef)^2)
print(as.data.frame(mse))

```


