#'@title simulate design matrix with a positive 
#'exchangeable correlation from 
#'a multivariate normal distribution
#'@param n positive integer for sample size
#'@param p positive integer for number of predictors
#'@param pos_corr positive value in (0,1) for positive correlation coefficient
#'@param sd positive numeric scalar for the standard deviation 
#'@return a covariate matrix of dimension n by p
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#' X <- simul_covar_mat(n = 1000, p = 5, pos_corr = 0.5, sd = 2)
#' print(cor(X))
#' print(cov(X))
#'@export
simul_covar_mat <- function(n, p, pos_corr=0.5, sd = 1){
  if(pos_corr < 1){
    X <- matrix(rnorm(n*p,sd=sd*sqrt(1-pos_corr)), nrow=n,ncol=p)
    if(pos_corr > 0 ){
      Z <- rnorm(n,sd=sd*sqrt(pos_corr))
      X <- X + Z
    }
  } else {
    X <- matrix(rnorm(n,sd = sd),nrow = n, ncol = p)
  }
  return(X)
}


#'@title simulate data for a linear regression
#'@param X numeric matrix of dimension
#'@param regcoef numeric vector for regression coefficients
#'@param sd_eps positive numeric scalar for the random noise standard deviation 
#'@param R_sq positive value in (0,1) for coefficient of 
#'determination (R-squared) if sd_eps = NULL
#'@param n positive integer for sample size if X is NULL
#'@param p positive integer for number of covariates if X is NULL
#'@param pos_corr positive value in (0,1) for positive correlation between 
#'covariates if X is NULL
#'@param sd_X positive numeric scalar for the standard deviation of 
#'covariates if X is NULL
#'@return list object of 
#'\describe{
#'\item{y}{vector of n response variables}
#'\item{X}{covariate matrix}
#'\item{regcoef}{vector of p regression coefficients}
#'\item{intercept}{numeric scalar}
#'\item{mu}{vector of E(y | X)}
#'\item{sd_eps}{noise standard deviation}
#'\item{R_sq}{R squared}
#'}
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#'regdat <- simul_linear_reg(n = 1000, p = 5, pos_corr = 0.6)
#'res <- with(regdat, lm(y~X))
#'print(cbind(true_coef = c(regdat$intercept,regdat$regcoef), 
#'            est_coef = coef(res)))
#'@export
simul_linear_reg <- function(regcoef = NULL, 
                             X = NULL, 
                             intercept = 0,
                             sd_eps = NULL,
                             R_sq = 0.9,
                             n = 1000, 
                             p = length(regcoef), 
                             pos_corr=0.5, 
                             sd_X = 1){
  if(is.null(regcoef)){
    regcoef <- rep(1, length = p)
  }
  p = length(regcoef)
  if(is.null(X)){
    X <- simul_covar_mat(n, p, pos_corr, sd_X)
  }
  mu <- intercept + X%*%regcoef
  if(is.null(sd_eps)){
    sd_eps = sqrt((1.0 - R_sq)/R_sq)*sd(mu)
  }
  n <- nrow(X)
  y = mu + rnorm(n,sd = sd_eps)
  return(list(y = y, X = X, regcoef = regcoef, 
              intercept = intercept, mu = mu, 
              sd_eps = sd_eps, R_sq = R_sq, n = n, p = p))
}

#'@title simulate Gaussian mixture regression model
#'@param n positive integer for sample size
#'@param p positive integer for number of covariates
#'@param K positive integer for number of components
#'@param pos_corr numeric scalar in (0,1) for exchangable correlation between covariates
#'@param min_effect positive numeric scalar for the miminum effect size of the regression coefficients
#'@param max_effect positive numeric scalar for the maximum effect size of the regression coefficients
#'@param prec_prob positive numeric scalar control the weight probabilty parameters 
#'@param R_sq numeric scalar in (0,1) for R squared
#'@author Jian Kang <jiankang@umich.edu>
#'@return list object of 
#'\describe{
#'\item{y}{vector of n response variables}
#'\item{X}{covariate matrix}
#'\item{alpha}{p by K matrix for regression coefficients}
#'\item{prob}{K by 1 vector for weight probabilit}
#'\item{sigma_sq}{positive numeric scalar for the noise variance}
#'\item{delta}{n by 1 vector for group labels}
#'\item{R_sq}{R squared}
#'}
#'
#'@examples
#'library(flexmix)
#'dat <- simul_mixture_reg(n = 1000, K = 3)
#'fit <- flexmix(dat$y~0+dat$X,k = 3)
#'paras <- parameters(fit)[1:ncol(dat$X),]
#'print(paras) 
#'@export
simul_mixture_reg <- function(n = 1000, 
                              p = 5, 
                              K = 3, 
                              pos_corr = 0.5, 
                              min_effect = 1,
                              max_effect = 3, 
                              prec_prob = 50, 
                              R_sq = 0.9){
  X <- simul_covar_mat(n=n, p=p, pos_corr=pos_corr)
  alpha <- matrix(sample(c(-1,0,1),K*p,replace=TRUE)*runif(K*p,min_effect,max_effect),nrow=p,ncol=K)
  prob <- rgamma(K,shape=prec_prob,rate=1)
  prob <- prob/sum(prob)
  delta <- sample(K,n,replace=TRUE,prob=prob)
  mu <- rep(NA, length=n)
  for(k in 1:K){
    idx <- which(delta==k)
    mu[idx] = X[idx,]%*%alpha[,k]
  }
  sigma_sq = var(mu)*(1 - R_sq)/R_sq
  y = mu + rnorm(n,sd=sqrt(sigma_sq))
  return(list(y=y,X=X,alpha=alpha, prob=prob,sigma_sq=sigma_sq, delta = delta,R_sq=R_sq))
}


#'@title find the column index permutation to match two matrix parameters
#'@param target_paras matrix of parameters
#'@param paras matrix of parameters 
#'@return column indices idx_match such that target_paras is matched with paras[,idx_match]
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#'library(flexmix)
#'set.seed(2024)
#'dat <- simul_mixture_reg(n = 1000, K = 3)
#'fit <- flexmix(dat$y~0+dat$X,k = 3)
#'alpha_est <- parameters(fit)[1:ncol(dat$X),]
#'print(alpha_est) 
#'idx <- match_paras(dat$alpha,alpha_est)
#'print(dat$alpha)
#'print(alpha_est[,idx])
#'print(idx)
#'@export
match_paras <- function(target_paras,paras){
  idx_match = rep(NA,length=ncol(target_paras))
  for(i in 1:ncol(target_paras)){
    idx_match[i] <- which.min(apply((paras - target_paras[,i])^2,2,sum))
  }
  return(idx_match)
}

#'@title Fit ridge regression by different methods
#'@param y numeric vector for n observations of the response variable
#'@param X covariate matrix of dimension n by p
#'@param intercept logical variable indicating whether include the intercept 
#'in the model. The default value is TRUE
#'@param lambda tuning parameter 
#'@param method character indicating one of the following methods 
#'c("solve","solve-large-p","gd","gd-exact","gd-back","gd-BB"). 
#'"solve": solving linear system of dimensions p by p
#'"solve-large-p": solving linear system of dimensions n by n
#'"gd": Gradient descent with a fixed step size
#'"gd-exact": Gradient descent with exact line search 
#'"gd-back": Gradient descent with backtracking line search
#'"gd-BB": Gradient descent with Barzilai Borwein method to update the step size
#'"sgd": Stochastic gradient descent with a fixed step size
#'@return list objects with the following elements
#'\describe{
#'\item{intercept}{intercept estimate}
#'\item{regcoef}{a vector of p regression coefficient estimates}
#'\item{residuals}{a vector of n residuals}
#'\item{fitted.values}{the fitted mean values}
#'\item{gd}{list object including convergence information}
#'\item{elapsed}{computing time in seconds}
#'\item{loss}{value of loss function at the parameter estimates}
#'}
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#'library(glmnet)
#'set.seed(2023)
#'betacoef <- c(rep(c(1,-1),length=10),rep(0,length=200 - 10))
#'dat <- simul_linear_reg(regcoef = betacoef,intercept=2,
#' n = 10000,R_sq = 0.9)
#'lambda <- 0.05
#'method_list <- c("solve","solve-large-p","gd","gd-exact","gd-back","gd-BB","sgd")
#'res_list <- list()
#'loss_list <- list()
#'mse_list <- list()
#'elapsed_list <- list()
#'for(method in method_list){
#'   res_list[[method]] <- with(dat, linear_ridge_reg(y = y, X = X, lambda = lambda, method=method))
#'   loss_list[[method]] <- res_list[[method]]$loss
#'   mse_list[[method]] <- mean((res_list[[method]]$regcoef - dat$regcoef)^2)
#'   elapsed_list[[method]] <- as.numeric(res_list[[method]]$elapsed)
#'}
#'elapsed_list[["glmnet"]] <- proc.time()[3]
#'res_list[["glmnet"]] <- with(dat,glmnet(x=X, y=y, alpha=0.0, lambda = lambda))
#'elapsed_list[["glmnet"]] <- proc.time()[3] - elapsed_list[["glmnet"]]
#'loss_list[["glmnet"]] <- with(res_list[["glmnet"]],ridge_reg_loss(intercept=a0,regcoef=beta,lambda=lambda,dat$y,dat$X))
#'mse_list[["glmnet"]] <- mean((res_list[["glmnet"]]$beta - dat$regcoef)^2)
#'loss_list[["true"]] <- with(dat,ridge_reg_loss(intercept,regcoef,lambda=lambda,y,X))
#'mse_list[["true"]] <- NA
#'elapsed_list[["true"]] <- NA
#'print(cbind(loss = unlist(loss_list),
#'            mse = unlist(mse_list),
#'            elapsed = unlist(elapsed_list)))
#'@export
linear_ridge_reg <- function(y, X, 
                             intercept = TRUE, 
                             lambda = 0.001, 
                             method = c("solve","solve-large-p",
                                        "gd","gd-exact","gd-back","gd-BB",
                                        "sgd","sgd-decrease","sgd-batch"),
                             initial_regcoef = rep(0,length=ncol(X)),
                             initial_intercept = 0,
                             step_size = 1e-4,
                             max_iter = 100000L,
                             tol = 1e-4,
                             alpha = 0.25,
                             beta = 0.5,
                             diminishing_ratio = 0.1,
                             r = 0.55,
                             batch_size = 10,
                             verbose = Inf,
                             comp_loss = FALSE){
  elapsed <- proc.time()[3]
  if(intercept){
    X1 <- cbind(1, X)
    initial_reg1coef <- c(initial_intercept, initial_regcoef)
  } else{
    X1 <- X
    initial_reg1coef <- initial_regcoef
  }
  n = nrow(X1)
  p = ncol(X1)
  nlambda <- n*lambda
  method = method[1]
  gd = NULL
  loss_trace = NULL
  if(comp_loss==TRUE){
    loss_trace = rep(NA,length=max_iter)
  } 
  if(method == "solve"){
    #solve p by p matrix
    XtX <- crossprod(X1)
    Xty <- crossprod(X1,y)
    if(intercept){
      diag(XtX) <- diag(XtX) + c(0,rep(nlambda,length=p-1))
    } else{
      diag(XtX) <- diag(XtX) + nlambda
    }
    reg1coef <- solve(XtX,Xty)
  }
  if(method == "solve-large-p"){
    #solve n by n matrix
    if(intercept){
      XXt <- tcrossprod(X)
      diag(XXt) <- diag(XXt) + nlambda
      y1 = cbind(y,1)
      reg1coef <-  solve(XXt,y1)
      reg1coef <- crossprod(X,reg1coef)
      Xbeta <- X%*%reg1coef
      intercept <- mean(y - Xbeta[,1])/(1 - mean(Xbeta[,2]))
      regcoef <- reg1coef[,1] - intercept*reg1coef[,2]
      reg1coef <- c(intercept,regcoef)
    } else{
      XXt <- tcrossprod(X1)
      diag(XXt) <- diag(XXt) + nlambda
      reg1coef <-  solve(XXt,y)
      reg1coef <- crossprod(X1,reg1coef)
    }
  }
  if(method == "gd"){
    #gradient descent with fixed step size
    reg1coef <- initial_reg1coef
    iter = 0
    eps = 1
    Xbeta <- X1%*%reg1coef
    resid <- y - Xbeta
    
    while(iter < max_iter & eps > tol & eps < 1e10){
      Delta_beta <- crossprod(X1,resid)/n
      if(intercept){
        Delta_beta <- Delta_beta - c(0,lambda*reg1coef[-1])
      }
      else {
        Delta_beta <- Delta_beta - lambda*reg1coef
      }
      reg1coef <- reg1coef + step_size*Delta_beta
      Xbeta <- X1%*%reg1coef
      resid <- y - Xbeta
      #eps <- sum(step_size*abs(Delta_beta))
      eps <- mean(Delta_beta^2)
      iter = iter + 1
      if(iter%%verbose==0){
        if(intercept){
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
        cat(iter, loss, eps,"\n")
      }
      if(comp_loss){
        if(intercept){
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
      }
    }
    
    if(comp_loss){
      loss_trace = loss_trace[1:iter]
    }
    gd = list(iter=iter,eps=eps,convergence=as.numeric(eps < tol),
              step_size=step_size,
              loss_trace=loss_trace)
  }
  
  if(method == "gd-back"){
    #gradient descent with backtracking line search
    reg1coef = initial_reg1coef
    iter = 0
    eps = 1
    if(intercept){
      f <- function(x) return(ridge_reg_loss(x[1],x[-1],lambda,y,X))
    } else{
      f <- function(x) return(ridge_reg_loss(0,x,lambda,y,X))
    }
    while(iter < max_iter & eps > tol & eps < 1e10 & step_size!=0){
      Xbeta <- X1%*%reg1coef
      resid <- y - Xbeta
      Delta_beta <- crossprod(X1, resid)/n
      if(intercept){
        Delta_beta <- Delta_beta - c(0,lambda*reg1coef[-1])
      } else{
        Delta_beta <- Delta_beta - lambda*reg1coef
      }
      step_size <- backtracking_line_search(x=reg1coef, grad_x = -Delta_beta, f=f)
      reg1coef <- reg1coef + step_size*Delta_beta
      #eps <- sum(step_size*abs(Delta_beta))
      eps <- mean(Delta_beta^2)
      iter = iter + 1
      if(iter%%verbose==0){
        if(intercept)
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        else
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        cat(iter, loss, eps,"\n")
      }
      if(comp_loss){
        if(intercept){
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
      }
    }
    if(comp_loss){
      loss_trace = loss_trace[1:iter]
    }
    gd = list(iter=iter,eps=eps,convergence=as.numeric(eps < tol),
              step_size=step_size,
              loss_trace=loss_trace)
  }
  
  if(method == "gd-exact"){
    #gradient descent with exact line search
    reg1coef = initial_reg1coef
    iter = 0
    eps = 1
    
    while(iter < max_iter & eps > tol & eps < 1e10 & step_size!=0){
      Xbeta <- X1%*%reg1coef
      resid <- y - Xbeta
      Delta_beta <- crossprod(X1, resid)/n
      if(intercept){
        Delta_beta <- Delta_beta - c(0,lambda*reg1coef[-1])
      } else{
        Delta_beta <- Delta_beta - lambda*reg1coef
      }
      step_size <-  ridge_reg_exact_line_search(reg1coef, Delta_beta, lambda, X1, resid)
      
      reg1coef <- reg1coef + step_size*Delta_beta
      eps <- sum(step_size*abs(Delta_beta))
      iter = iter + 1
      if(iter%%verbose==0){
        if(intercept){
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
        cat(iter, loss, eps,"\n")
      }
      if(comp_loss){
        if(intercept){
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
      }
      
    }
    if(comp_loss){
      loss_trace = loss_trace[1:iter]
    }
    gd = list(iter=iter,eps=eps,convergence=as.numeric(eps < tol),
              step_size=step_size,
              loss_trace=loss_trace)
  }
  
  if(method == "gd-BB"){
    #gradient descent with BB methods
    reg1coef = initial_reg1coef
    iter = 0
    eps = 1
    
    Xbeta <- X1%*%reg1coef
    resid <- y - Xbeta
    
    Delta_beta <- crossprod(X1, resid)/n
    if(intercept){
      Delta_beta <- Delta_beta - c(0,lambda*reg1coef[-1])
    } else{
      Delta_beta <- Delta_beta - lambda*reg1coef
    }
    prev_reg1coef = reg1coef
    reg1coef <- reg1coef + step_size*Delta_beta
    
    while(iter < max_iter & eps > tol & eps < 1e10 & step_size!=0){
      Xbeta <- X1%*%reg1coef
      resid <- y - Xbeta
      
      prev_Delta_beta <- Delta_beta
      
      Delta_beta <- crossprod(X1, resid)/n
      if(intercept){
        Delta_beta <- Delta_beta - c(0,lambda*reg1coef[-1])
      } else{
        Delta_beta <- Delta_beta - lambda*reg1coef
      }
      step_size <- Barzilai_Borwein_method(reg1coef, prev_reg1coef,-Delta_beta, -prev_Delta_beta)
      
      prev_reg1coef <- reg1coef
      reg1coef <- reg1coef + step_size*Delta_beta
      eps <- sum(step_size*abs(Delta_beta))
      iter = iter + 1
      if(iter%%verbose==0){
        if(intercept){
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
        cat(iter, loss, eps,"\n")
      }
      if(comp_loss){
        if(intercept){
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
      }
      
    }
    if(comp_loss){
      loss_trace = loss_trace[1:iter]
    }
    gd = list(iter=iter,eps=eps,convergence=as.numeric(eps < tol),
              step_size=step_size,
              loss_trace=loss_trace)
  }
  
  
  if(method == "sgd"){
    #stochastic gradient descent with fixed step size
    reg1coef <- initial_reg1coef
    iter = 0
    eps = 1
    Xbeta <- X1%*%reg1coef
    resid <- y - Xbeta
    n <- length(y)
    while(iter < max_iter & eps < 1e10){
      i = sample(1:n,1)
      Delta_beta <- X1[i,]*(y[i] - sum(X1[i,]*reg1coef))
      
      if(intercept){
        Delta_beta <- Delta_beta - c(0,lambda*reg1coef[-1])
      }
      else {
        Delta_beta <- Delta_beta - lambda*reg1coef
      }
      reg1coef <- reg1coef + step_size*Delta_beta
      iter = iter + 1
      if(iter%%verbose==0){
        if(intercept){
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
        cat(iter, loss, eps,"\n")
      }
      if(comp_loss){
        Xbeta <- X1%*%reg1coef
        resid <- y - Xbeta
        if(intercept){
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
      }
    }
    
    if(comp_loss){
      loss_trace = loss_trace[1:iter]
    }
    gd = list(iter=iter,
              step_size=step_size,
              loss_trace=loss_trace)
  }
  
  if(method == "sgd-decrease"){
    #stochastic gradient descent with diminish
    b_0 <- max_iter/((1.0/diminishing_ratio)^(1.0/r) - 1.0)
    a_0 <- step_size*b_0^r
    reg1coef <- initial_reg1coef
    iter = 0
    eps = 1
    Xbeta <- X1%*%reg1coef
    resid <- y - Xbeta
    n <- length(y)
    while(iter < max_iter & eps < 1e10){
      i = sample(1:n,1)
      Delta_beta <- X1[i,]*(y[i] - sum(X1[i,]*reg1coef))
      
      if(intercept){
        Delta_beta <- Delta_beta - c(0,lambda*reg1coef[-1])
      }
      else {
        Delta_beta <- Delta_beta - lambda*reg1coef
      }
      step_size = a_0*(b_0 + iter)^(-r)
      reg1coef <- reg1coef + step_size*Delta_beta
      iter = iter + 1
      if(iter%%verbose==0){
        if(intercept){
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
        cat(iter, loss, eps,"\n")
      }
      if(comp_loss){
        Xbeta <- X1%*%reg1coef
        resid <- y - Xbeta
        if(intercept){
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
      }
    }
    
    gd = list(iter=iter,
              step_size=step_size,
              loss_trace=loss_trace)
  }
  
  if(method == "sgd-batch"){
    #stochastic gradient descent with diminish
    b_0 <- max_iter/((1.0/diminishing_ratio)^(1.0/r) - 1.0)
    a_0 <- step_size*b_0^r
    reg1coef <- initial_reg1coef
    iter = 0
    eps = 1
    Xbeta <- X1%*%reg1coef
    resid <- y - Xbeta
    n <- length(y)
    while(iter < max_iter & eps < 1e10){
      I_b = sample(1:n,batch_size)
      Delta_beta <- crossprod(X1[I_b,],(y[I_b] - X1[I_b,]%*%reg1coef))/batch_size
      
      if(intercept){
        Delta_beta <- Delta_beta - c(0,lambda*reg1coef[-1])
      }
      else {
        Delta_beta <- Delta_beta - lambda*reg1coef
      }
      step_size = a_0*(b_0 + iter)^(-r)
      reg1coef <- reg1coef + step_size*Delta_beta
      iter = iter + 1
      if(iter%%verbose==0){
        if(intercept){
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
        cat(iter, loss, eps,"\n")
      }
      if(comp_loss){
        Xbeta <- X1%*%reg1coef
        resid <- y - Xbeta
        if(intercept){
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2))
        } else{
          loss_trace[iter] <- 0.5*(mean(resid^2) + lambda*sum(reg1coef^2))
        }
      }
    }
    
    gd = list(iter=iter,
              step_size=step_size,
              loss_trace=loss_trace)
  }
  
  fitted.values <- X1%*%reg1coef
  resid <- y - fitted.values
  if(intercept){
    res <- list(intercept = reg1coef[1],
                regcoef = reg1coef[-1],
                residuals = resid,
                fitted.values = fitted.values,
                lambda = lambda,
                method = method,
                loss = 0.5*(mean(resid^2) + lambda*sum(reg1coef[-1]^2)),
                gd = gd,
                elapsed = proc.time()[3] - elapsed)
  } else{
    res <- list(intercept = 0,
                regcoef = reg1coef,
                residuals = resid,
                fitted.values = fitted.values,
                lambda = lambda,
                method = method,
                loss = 0.5*(mean(resid^2) + lambda*sum(reg1coef^2)),
                gd = gd,
                elapsed = proc.time()[3] - elapsed)
  }
  return(res)
}

#'@title Evaluate loss function for ridge regression
#'@param intercept numeric scalar
#'@param regcoef numeric vector of length p
#'@param lambda numeric scalar tuning parameter
#'@param y numeric vector for n observations of the response variable
#'@param X covariate matrix of dimension n by p
#'@return value of the loss function ||y - X*beta||^2_2/n + lambda*||X beta||^2_2/2
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#'library(glmnet)
#'set.seed(2023)
#'regdat <- simul_linear_reg(regcoef = c(rep(c(1,-1),length=10),rep(0,length=990)),
#'intercept = 2, 
#'n = 100, R_sq = 0.9)
#'lambda <- 0.04
#'res_glmnet <- with(regdat,glmnet(x=X,y=y,lambda=lambda,alpha=0.0))
#'res_solve <- with(regdat, linear_ridge_reg(y = y, X = X, lambda = lambda,method="solve"))
#'res_solve_p <- with(regdat, linear_ridge_reg(y = y, X = X, lambda = lambda,method="solve-large-p"))
#'loss <- list()
#'loss$glmnet <- with(res_glmnet,ridge_reg_loss(intercept=a0,regcoef=beta,lambda=lambda,regdat$y,regdat$X))
#'loss$solve <- with(res_solve,ridge_reg_loss(intercept,regcoef,lambda=lambda,regdat$y,regdat$X))
#'loss$solve_p <- with(res_solve_p,ridge_reg_loss(intercept,regcoef,lambda=lambda,regdat$y,regdat$X))
#'loss$true <- with(regdat,ridge_reg_loss(intercept,regcoef,lambda=lambda,y,X))
#'print(as.data.frame(loss))
#'mse <- list()
#'mse$glmnet <- mean((res_glmnet$beta - regdat$regcoef)^2)
#'mse$solve <- mean((res_solve$regcoef - regdat$regcoef)^2)
#'mse$solve_p <- mean((res_solve_p$regcoef - regdat$regcoef)^2)
#'print(as.data.frame(mse))
#'@export
ridge_reg_loss <- function(intercept, regcoef, lambda, y, X){
  resid <- y - X%*%regcoef - intercept
  loss <- 0.5*(mean(resid*resid) + lambda*sum(regcoef*regcoef))
  return(loss)
}

#'@title Backtracking line search for gradient decent
#'@param x numeric vector for x in dom f
#'@param grad_x numeric vector for gradient f
#'@param f R function for the objective function 
#'@param alpha numeric scalar between (0, 0.5)
#'@param beta numeric scalar between (0, 1)
#'@return numeric value for the step_size
#'@author Jian Kang <jiankang@umich.edu>
#'@examples
#'set.seed(2023)
#'regdat <- simul_linear_reg(regcoef = c(rep(c(1,-1),length=10),rep(0,length=90)),
#'intercept = 2, 
#'n = 50, R_sq = 0.9)
#'lambda <- 0.04
#'n = length(regdat$y)
#'p = ncol(regdat$X)
#'b <- rnorm(p+1)
#'X1 <- cbind(1,regdat$X)
#'grad_b <- crossprod(X1, X1%*%b - regdat$y)/n + lambda*b
#'f <- function(x) return(ridge_reg_loss(x[1],x[-1],lambda,regdat$y,regdat$X))
#'step_size <-  backtracking_line_search(b,grad_b,f=f)
#'print(c(f(b),f(b-step_size*grad_b)))
#'@export
backtracking_line_search <- function(x, grad_x, f, alpha = 0.25, beta = 0.5, max_iter=1000L, ...){
  step_size = 1
  sum_grad_sq<- sum(grad_x*grad_x)
  fx <- f(x,...)
  iter = 0
  while(f(x - step_size*grad_x,...) > fx - alpha*step_size*sum_grad_sq & iter < max_iter){
    step_size <- step_size*beta
    iter = iter + 1
  }
  return(step_size)
}


#'@title Exact line search for ridge regression
#'@param beta regression coefficients
#'@param Delta_beta descent direction for regression coefficients
#'@param lambda tuning parameter
#'@param X design matrix
#'@param resid residuals
#'@return step_size
#'@author Jian Kang <jiankang@umich.edu>
#'@examples 
#'set.seed(2023)
#'regdat <- simul_linear_reg(regcoef = c(rep(c(1,-1),length=10),rep(0,length=90)),
#'intercept = 2, 
#'n = 50, R_sq = 0.9)
#'lambda <- 0.04
#'n = length(regdat$y)
#'p = ncol(regdat$X)
#'b <- rnorm(p+1)
#'X1 <- cbind(1,regdat$X)
#'resid <- regdat$y - X1%*%b
#'grad_b <- crossprod(X1,-resid)/n + lambda*b
#'f <- function(x) return(ridge_reg_loss(x[1],x[-1],lambda,regdat$y,regdat$X))
#'step_size <-  ridge_reg_exact_line_search(b, -grad_b, lambda, X1, resid)
#'print(c(f(b),f(b-step_size*grad_b)))
#'@export
ridge_reg_exact_line_search <- function(beta,Delta_beta,lambda,X,resid,step_size_0=0){
  X_Delta_beta <- X%*%Delta_beta
  step_size <- sum(resid*X_Delta_beta) + lambda*sum(beta*Delta_beta)
  if(step_size > step_size_0){
    step_size <- step_size/(sum(X_Delta_beta^2) + lambda*sum(Delta_beta^2))
  } else{
    step_size = step_size_0
  }
  return(step_size)
}

#'@title Barzilai Borwein method to determine the step size for gradient descent 
#'@param x  numeric vector for current parameters
#'@param prev_x numeric vector for previous parameters
#'@param grad_x  numeric vector for current gradient
#'@param prev_grad_x vector for previous gradient
#'@param step_size lower bound of the step size
#'@return numeric scalar for step size 
#'@author Jian Kang <jiankang@umich.edu>
#'@examples 
#'set.seed(2023)
#'regdat <- simul_linear_reg(regcoef = c(rep(c(1,-1),length=5),rep(0,length=5)),
#'intercept = 2, 
#'n = 500, R_sq = 0.9)
#'lambda <- 0.04
#'step_size <- 0.001
#'n = length(regdat$y)
#'p = ncol(regdat$X)
#'b <- rnorm(p+1)
#'X1 <- cbind(1,regdat$X)
#'resid <- regdat$y - X1%*%b
#'grad_b <- crossprod(X1,-resid)/n + lambda*b
#'prev_b <- b
#'b <- b - step_size*grad_b
#'resid <- regdat$y - X1%*%b
#'prev_grad_b <- grad_b
#'grad_b <- crossprod(X1,-resid)/n + lambda*b
#'f <- function(x) return(ridge_reg_loss(x[1],x[-1],lambda,regdat$y,regdat$X))
#'prev_step_size <- step_size
#'step_size <-  Barzilai_Borwein_method(b, prev_b, grad_b, prev_grad_b)
#'print(c(f(b+prev_step_size*prev_grad_b),f(b),f(b-step_size*grad_b)))
#'@export
Barzilai_Borwein_method <- function(x, prev_x, grad_x, prev_grad_x,step_size=1e-6){
  diff_x = x - prev_x
  diff_grad_x = grad_x - prev_grad_x
  step_size = max(sum(diff_x*diff_grad_x)/sum(diff_grad_x*diff_grad_x),step_size)
  #step_size = max(sum(diff_x^2)/sum(diff_x*diff_grad_x), step_size)
  return(step_size)
}