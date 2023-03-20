#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//'@useDynLib UmichBiostat815, .registration=TRUE

using namespace Rcpp;
using namespace arma;

//'@title Barzilai Borwein method to determine the step size for gradient descent 
//'@param x  numeric vector for current parameters
//'@param prev_x numeric vector for previous parameters
//'@param grad_x  numeric vector for current gradient
//'@param prev_grad_x vector for previous gradient
//'@param step_size lower bound of the step size
//'@return numeric scalar for step size 
//'@author Jian Kang <jiankang@umich.edu>
//'@examples 
//'set.seed(2023)
//'regdat <- simul_linear_reg(regcoef = c(rep(c(1,-1),length=5),rep(0,length=5)),
//'intercept = 2, 
//'n = 500, R_sq = 0.9)
//'lambda <- 0.04
//'step_size <- 0.001
//'n = length(regdat$y)
//'p = ncol(regdat$X)
//'b <- rnorm(p+1)
//'X1 <- cbind(1,regdat$X)
//'resid <- regdat$y - X1%*%b
//'grad_b <- crossprod(X1,-resid)/n + lambda*b
//'prev_b <- b
//'b <- b - step_size*grad_b
//'resid <- regdat$y - X1%*%b
//'prev_grad_b <- grad_b
//'grad_b <- crossprod(X1,-resid)/n + lambda*b
//'f <- function(x) return(ridge_reg_loss(x[1],x[-1],lambda,regdat$y,regdat$X))
//'prev_step_size <- step_size
//'step_size <-  Barzilai_Borwein_method_cpp(b, prev_b, grad_b, prev_grad_b)
//'print(c(f(b+prev_step_size*prev_grad_b),f(b),f(b-step_size*grad_b)))
//'@export
// [[Rcpp::export]]
double Barzilai_Borwein_method_cpp(arma::vec& x, 
                                   arma::vec& prev_x, 
                                   arma::vec& grad_x, 
                                   arma::vec& prev_grad_x,
                                   double step_size_lb =1e-6){
  arma::vec diff_x = x - prev_x;
  arma::vec diff_grad_x = grad_x - prev_grad_x;
  double step_size = arma::sum(diff_x%diff_x)/arma::sum(diff_x%diff_grad_x);
  if(step_size < 0){
    step_size = step_size_lb;
  } 
  return step_size;
}


//'@title Backtraking line search in Rcpp
//'@param paras current value of the parameter
//'@param grad_paras current value of the gradient parameter
//'@param obj_fun objective function to be minimized
//'@param dat list object for additional information for the objective function
//'@param a numeric scalar between (0, 0.5)
//'@param b numeric scalar between (0, 1)
//'@return numeric value for the step_size
//'@author Jian Kang <jiankang@umich.edu>
//'@export
// [[Rcpp::export]]
double backtracking_line_search_cpp(arma::vec& paras, 
                                     arma::vec& grad_paras, 
                                     Rcpp::Function obj_fun, 
                                     Rcpp::List& dat,
                                     double a = 0.25, 
                                     double b = 0.5, 
                                     int max_iter = 1000L){
   double step_size = 1;
   double sum_grad_sq =  arma::sum(grad_paras%grad_paras);
   double f_paras_plus = Rcpp::as<double>(obj_fun(paras - step_size*grad_paras,dat));
   double f_paras = Rcpp::as<double>(obj_fun(paras,dat));
   int iter= 0;
   while(f_paras_plus > f_paras - a*step_size*sum_grad_sq && iter < max_iter){
     step_size *= b;
     f_paras_plus = Rcpp::as<double>(obj_fun(paras - step_size*grad_paras,dat));
     iter++;
   }
   return(step_size);
}



//'@title Gradient descent algorithm with different step size methods
//'@param paras numeric vector
//'@param obj_fun objective function to be minimized
//'@param grad_fun gradient function
//'@param dat list object to pass other data to the objective function 
//'and the gradient of the objective functions
//'@param step_size_method indicate the methods to determine the step size
//'@param max_iter positive integer for maximum number of iterations
//'@param step_size numeric scalar for the step size
//'@param tol tolerance for the stopping criterion
//'@author Jian Kang <jiankang@umich.edu>
//'@return a list of objects
//'@examples 
//'regdat <- simul_linear_reg(n = 1000, p = 5, pos_corr = 0.6)
//'regdat$lambda = 0.001
//'res <- gradient_descent_cpp(paras=rep(0,length=ncol(regdat$X)),
//'ridge_reg_loss_cpp,grad_ridge_reg_loss_cpp,regdat)
//'@export
// [[Rcpp::export]]
Rcpp::List gradient_descent_cpp(arma::vec paras, 
                      Rcpp::Function obj_fun, 
                      Rcpp::Function grad_fun, 
                      Rcpp::List& dat,
                      Rcpp::CharacterVector step_size_method = Rcpp::CharacterVector::create("fixed","back","BB"),
                      int max_iter = 100000L,
                      double step_size = 1e-4,
                      double tol = 1e-4,
                      int verbose = 100L)
{
    arma::wall_clock timer;
    timer.tic();
    step_size_method = step_size_method[0];
    arma::vec value = Rcpp::as<arma::vec>(Rcpp::wrap(obj_fun(paras, dat)));
    arma::vec initial_paras = paras;
    arma::vec initial_value = value;
    int convergence = 0;
    double eps = 1;
    int iter = 0;
    arma::vec prev_grad_value; 
    arma::vec prev_paras; 
    arma::vec grad_value;
    NumericVector grad_value_0;
    if(step_size_method[0]=="BB"){
      grad_value_0 = grad_fun(paras, dat);
      grad_value = arma::vec(grad_value_0.begin(), grad_value_0.length(),false);
      prev_paras = paras;
      paras -= step_size*grad_value;
    }
    while(iter<max_iter && eps > tol && eps < 1e5){
      grad_value_0 = grad_fun(paras, dat);
      if(step_size_method[0]=="BB"){
        prev_grad_value = grad_value;
      }
      grad_value = arma::vec(grad_value_0.begin(), grad_value_0.length(),false);
      eps = mean(grad_value%grad_value);
      if(step_size_method[0]=="back"){
        step_size = backtracking_line_search_cpp(paras,grad_value,obj_fun,dat);
      } 
      if(step_size_method[0]=="BB"){
        step_size = Barzilai_Borwein_method_cpp(paras,prev_paras,grad_value,prev_grad_value);
        prev_paras = paras;
      }
      paras -= step_size*grad_value;
      iter++;
      if(iter%verbose==0){
        std::cout << iter << ": " << eps << std::endl;
      }
    }
    value = Rcpp::as<arma::vec>(Rcpp::wrap(obj_fun(paras, dat)));
    if(iter<max_iter && eps < tol){
      convergence = 1;
    }
    double elapsed = timer.toc();
    Rcpp::List res = Rcpp::List::create(Rcpp::Named("paras") = paras,
                            Rcpp::Named("value") = value,
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("convergence") = convergence,
                            Rcpp::Named("step_size_method") = step_size_method,
                            Rcpp::Named("elapsed") = elapsed);

    return res ;
}

//'@title Stochastic Gradient descent algorithm with different step size methods
//'@param paras numeric vector
//'@param obj_fun objective function to be minimized
//'@param grad_fun gradient function
//'@param dat list object to pass other data to the objective function 
//'and the gradient of the objective functions
//'@param step_size_method indicate the methods to determine the step size
//'@param n total numbers of sub functions in the objective function
//'@param max_iter positive integer for maximum number of iterations
//'@param weight number between 0 and 1 to indicate the rate to smooth the gradient over iterations
//'@param step_size numeric scalar for the initial step size
//'@param diminishing_ratio the rate of decreasing the step size
//'@param tol tolerance for the stopping criterion
//'@author Jian Kang <jiankang@umich.edu>
//'@return a list of objects
//'@examples 
//'regdat <- simul_linear_reg(n = 1000, p = 5, pos_corr = 0.6)
//'regdat$lambda = 0.001
//'res <- sgd_cpp(paras=rep(0,length=ncol(regdat$X)),ridge_reg_loss_cpp,
//'grad_ridge_reg_loss_cpp,regdat,
//'step_size_method = "decreasing")
//'
//'@export
// [[Rcpp::export]]
Rcpp::List sgd_cpp(arma::vec paras, 
                   Rcpp::Function obj_fun, 
                   Rcpp::Function grad_fun, 
                   Rcpp::List& dat,
                   Rcpp::CharacterVector step_size_method = Rcpp::CharacterVector::create("decreasing","BB","fixed"),
                   int n = -1,
                   int max_iter = 10000L,
                   int mini_batch_size = 1L,
                   int update_freq = 100L,
                   double weight = 1.0,
                   double step_size = 1e-3,
                   double diminishing_ratio = 0.1,
                   double r = 0.55,
                   double tol = 1e-4,
                   int burnin = 5000L,
                   int verbose = 1000L)
{
  arma::wall_clock timer;
  timer.tic();
  step_size_method = step_size_method[0];
  arma::vec value = Rcpp::as<arma::vec>(Rcpp::wrap(obj_fun(paras, dat)));
  arma::vec initial_paras = paras;
  arma::vec initial_value = value;
  int convergence = 0;
  double eps = 1;
  double step_size_0 = step_size;
  int iter = 0;
  arma::vec current_grad_value;
  arma::vec prev_grad_value; 
  arma::vec mean_paras(paras.n_elem,arma::fill::zeros);
  arma::vec prev_paras; 
  arma::vec grad_value;
  arma::vec temp_grad_value;
  NumericVector grad_value_0;
  NumericVector temp_grad_value_0;
  if(n==-1){
    n = Rcpp::as<int>(dat["n"]);
    std::cout << n << std::endl;
  }
  Rcpp::IntegerVector idx;
  if(step_size_method[0]=="BB"){
     grad_value_0 = grad_fun(paras, dat);
     grad_value = arma::vec(grad_value_0.begin(), grad_value_0.length(),false);
     prev_paras = paras;
     paras -= step_size*grad_value;
     prev_grad_value = grad_value;
     current_grad_value.zeros(paras.n_elem);
  }
  
  double b_0;
  double a_0;
  if(step_size_method[0]=="decreasing"){
    b_0 = (max_iter*1.0)/(pow(1.0/diminishing_ratio,1.0/r) - 1.0);
    a_0 = step_size*pow(b_0,r);
  }
  int k = 0;
  while(iter<max_iter && eps > tol && eps < 1e5){
    idx = Rcpp::sample(n,mini_batch_size) - 1;
    grad_value_0 = grad_fun(paras, dat, idx);
    grad_value = arma::vec(grad_value_0.begin(), grad_value_0.length(),false);
    //std::cout << grad_value << std::endl;
    
    if(step_size_method[0]=="BB"){
      current_grad_value *= (1.0 - weight);
      current_grad_value += weight*grad_value;
    }
    if(step_size_method[0]=="BB" && iter%update_freq==0){
      eps = mean(current_grad_value%current_grad_value);
      if(eps>tol){
        step_size = Barzilai_Borwein_method_cpp(paras,
                                              prev_paras,
                                              current_grad_value,
                                              prev_grad_value);
        step_size /= update_freq;
        prev_paras = paras;
        prev_grad_value = current_grad_value;
        current_grad_value.zeros();
      }
    }
    paras -= step_size*grad_value;
    if(step_size_method[0]=="decreasing"){
      step_size =  a_0*pow(b_0 + iter, -r);
    } 
    if(iter>=burnin){
      mean_paras += paras;
      k++;
    }
    
    if(iter%verbose==0){
      temp_grad_value_0 = grad_fun(paras, dat);
      temp_grad_value = arma::vec(temp_grad_value_0.begin(), temp_grad_value_0.length(),false);
      eps = mean(temp_grad_value%temp_grad_value);
      std::cout << iter << ": " << eps << std::endl;
    }
    iter++;
  }
  value = Rcpp::as<arma::vec>(Rcpp::wrap(obj_fun(paras, dat)));
  grad_value_0 = grad_fun(paras, dat);
  grad_value = arma::vec(grad_value_0.begin(), grad_value_0.length(),false);
  eps = mean(grad_value%grad_value);
  if(iter<max_iter && eps < tol){
    convergence = 1;
  }
  if(k>0){
    mean_paras /= k;
  } else{
    mean_paras = paras;
  }
  double elapsed = timer.toc();
  Rcpp::List res = Rcpp::List::create(Rcpp::Named("paras") = mean_paras,
                                      Rcpp::Named("value") = value,
                                      Rcpp::Named("iter") = iter,
                                      Rcpp::Named("convergence") = convergence,
                                      Rcpp::Named("step_size_method") = step_size_method,
                                      Rcpp::Named("update_freq") = update_freq,
                                      Rcpp::Named("weight") = weight,
                                      Rcpp::Named("step_size") = step_size,
                                      Rcpp::Named("initial_step_size") = step_size_0,
                                      Rcpp::Named("elapsed") = elapsed) ;
  
  return res ;
}


//'@title Compute loss of the ridge regression in Rcpp
//'@param beta numeric vector for regression coefficients
//'@param dat list object containing covariate matrix, response variable
//'@param idx_R vector of positive integers for sample indices. 
//'@return value of the loss function
//'@author Jian Kang <jiankang@umich.edu>
//'@export
// [[Rcpp::export(rng = false)]]
double ridge_reg_loss_cpp(arma::vec& beta, Rcpp::List& dat,
                          Rcpp::IntegerVector idx_R = Rcpp::IntegerVector::create()){
  double loss = 0.0;
  //arma::mat X = Rcpp::as<arma::mat>(dat["X"]); 
  //arma::vec y = Rcpp::as<arma::mat>(dat["y"]);
  NumericVector y_temp = dat["y"];
  NumericMatrix X_temp = dat["X"];
  arma::mat X(X_temp.begin(), X_temp.nrow(), X_temp.ncol(), false);
  arma::vec y(y_temp.begin(), y_temp.length(),false);
  double lambda = Rcpp::as<double>(dat["lambda"]);
  arma::vec neg_resid;
  arma::uvec idx_arma;
  if(idx_R.length()==0){
    neg_resid = X*beta;
    neg_resid -= y;
  } else{
    idx_arma = Rcpp::as<arma::uvec>(idx_R);
    neg_resid = X.rows(idx_arma)*beta;
    neg_resid -= y.rows(idx_arma);
  }
  loss += arma::mean(neg_resid%neg_resid);
  loss += lambda*arma::sum(beta%beta);
  loss *= 0.5;
  return loss;
}



//'@title compute the gradient of the ridge regression loss in Rcpp
//'@param beta regression coefficients
//'@param dat list object containing covariate matrix, response variable
//'@param idx_R vector of positive integers for sample indices. 
//'@return a vector of gradient of the loss function
//'@author Jian Kang <jiankang@umich.edu>
//'@export
 // [[Rcpp::export(rng = false)]]
 arma::vec grad_ridge_reg_loss_cpp(arma::vec& beta, Rcpp::List& dat, 
                                   Rcpp::IntegerVector idx_R = Rcpp::IntegerVector::create()){
   //arma::mat X = Rcpp::as<arma::mat>(dat["X"]); 
   //arma::vec y = Rcpp::as<arma::mat>(dat["y"]);
   NumericVector y_temp = dat["y"];
   NumericMatrix X_temp = dat["X"];
   arma::mat X(X_temp.begin(), X_temp.nrow(), X_temp.ncol(), false);
   arma::vec y(y_temp.begin(), y_temp.length(),false);
   double lambda = Rcpp::as<double>(dat["lambda"]);
   arma::vec neg_resid;
   arma::uvec idx_arma;
   if(idx_R.length()==0){
     neg_resid = X*beta;
     neg_resid -= y;
     neg_resid /= y.n_elem;
   } else{
     idx_arma = Rcpp::as<arma::uvec>(idx_R);
     neg_resid = X.rows(idx_arma)*beta;
     neg_resid -= y.rows(idx_arma);
     neg_resid /= idx_arma.n_elem;
   }
   arma::vec grad;
   if(idx_R.length()==0){
     grad = X.t()*neg_resid;
   } else{
     grad = X.rows(idx_arma).t()*neg_resid;
   }
   grad += lambda*beta;
   return grad;
 }



