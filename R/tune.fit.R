#' Using the \pkg{glmnet}, \pkg{ncvreg}, \pkg{msaenet}, \pkg{Coxnet} and \pkg{gcdnet} packages, fits a Generalized Linear
#' Model or Cox Proportional Hazards Model using various methods for choosing
#' the regularization parameter \eqn{\lambda}
#'
#' This function fits a generalized linear model or a Cox proportional hazards
#' model via penalized maximum likelihood, with available penalties as
#' indicated in the \pkg{glmnet}, \pkg{ncvreg}, \pkg{msaenet}, \pkg{Coxnet} and \pkg{gcdnet} packages. Instead of
#' providing the whole regularization solution path, the function returns the
#' solution at a unique value of \eqn{\lambda}, the one optimizing the
#' criterion specified in \code{tune}.
#'
#' @importFrom ncvreg ncvsurv
#' @importFrom ncvreg cv.ncvsurv
#' @importFrom gcdnet cv.gcdnet
#' @importFrom glmnet cv.glmnet
#' @importFrom msaenet aenet
#' @import doParallel
#' @export
#' @param x The design matrix, of dimensions n * p, without an intercept. Each
#' row is an observation vector.
#' @param y The response vector of dimension n * 1. Quantitative for
#' \code{family='gaussian'}, non-negative counts for \code{family='poisson'},
#' binary (0-1) for \code{family='binomial'}. For \code{family='cox'}, \code{y}
#' should be an object of class \code{Surv}, as provided by the function
#' \code{Surv()} in the package \pkg{survival}.
#' @param family Response type (see above).
#' @param penalty The penalty to be applied in the regularized likelihood
#' subproblems. 'SCAD' (the default), 'MCP', or 'lasso' are provided.
#' @param concavity.parameter The tuning parameter used to adjust the concavity
#' of the SCAD/MCP penalty. Default is 3.7 for SCAD and 3 for MCP.
#' @param tune Method for selecting the regularization parameter along the
#' solution path of the penalized likelihood problem. Options to provide a
#' final model include \code{tune='cv'}, \code{tune='aic'}, \code{tune='bic'},
#' and \code{tune='ebic'}. See references at the end for details.
#' @param nfolds Number of folds used in cross-validation. The default is 10.
#' @param type.measure Loss to use for cross-validation. Currently five
#' options, not all available for all models. The default is
#' \code{type.measure='deviance'}, which uses squared-error for gaussian models
#' (also equivalent to \code{type.measure='mse'} in this case), deviance for
#' logistic and poisson regression, and partial-likelihood for the Cox model.
#' Both \code{type.measure='class'} and \code{type.measure='auc'} apply only to
#' logistic regression and give misclassification error and area under the ROC
#' curve, respectively. \code{type.measure='mse'} or \code{type.measure='mae'}
#' (mean absolute error) can be used by all models except the \code{'cox'};
#' they measure the deviation from the fitted mean to the response. For
#' \code{penalty='SCAD'}, \code{penalty='MCP'}, \code{penalty='aenet'} and \code{penalty='msaenet'}, only
#' \code{type.measure='deviance'} is available.
#' @param gamma.ebic Specifies the parameter in the Extended BIC criterion
#' penalizing the size of the corresponding model space. The default is
#' \code{gamma.ebic=1}. See references at the end for details.
#' @param parallel Specifies whether to conduct parallel computing
#' @param seed An optimal argument for setting the seed to ensure reproducibility
#' @return Returns an object with \item{ix}{ The vector of indices of the
#' nonzero coefficients selected by the maximum penalized likelihood procedure
#' with \code{tune} as the method for choosing the regularization parameter.  }
#' \item{a0}{The intercept of the final model selected by \code{tune}.  }
#' \item{beta}{The vector of coefficients of the final model selected by
#' \code{tune}.  }
#' \item{fit}{The fitted penalized regression object.}
#' \item{lambda}{The corresponding lambda in the final model.}
#' @author Jianqing Fan, Yang Feng, Diego Franco Saldana, Richard Samworth, Arce Domingo-Relloso and
#' Yichao Wu
#' @references Jerome Friedman and Trevor Hastie and Rob Tibshirani (2010)
#' Regularization Paths for Generalized Linear Models Via Coordinate Descent.
#' \emph{Journal of Statistical Software}, \bold{33}(1), 1-22.
#'
#' Noah Simon and Jerome Friedman and Trevor Hastie and Rob Tibshirani (2011)
#' Regularization Paths for Cox's Proportional Hazards Model Via Coordinate
#' Descent. \emph{Journal of Statistical Software}, \bold{39}(5), 1-13.
#'
#' Patrick Breheny and Jian Huang (2011) Coordiante Descent Algorithms for
#' Nonconvex Penalized Regression, with Applications to Biological Feature
#' Selection. \emph{The Annals of Applied Statistics}, \bold{5}, 232-253.
#'
#' Hirotogu Akaike (1973) Information Theory and an Extension of the Maximum
#' Likelihood Principle. In \emph{Proceedings of the 2nd International
#' Symposium on Information Theory}, BN Petrov and F Csaki (eds.), 267-281.
#'
#' Gideon Schwarz (1978) Estimating the Dimension of a Model. \emph{The Annals
#' of Statistics}, \bold{6}, 461-464.
#'
#' Jiahua Chen and Zehua Chen (2008) Extended Bayesian Information Criteria for
#' Model Selection with Large Model Spaces. \emph{Biometrika}, \bold{95},
#' 759-771.
#' @keywords models
#' @examples
#'
#'
#' set.seed(0)
#' data("leukemia.train", package = "SIS")
#' y.train <- leukemia.train[, dim(leukemia.train)[2]]
#' x.train <- as.matrix(leukemia.train[, -dim(leukemia.train)[2]])
#' x.train <- standardize(x.train)
#' model <- tune.fit(x.train[, 1:3500], y.train, family = "binomial", tune = "bic")
#' model$ix
#' model$a0
#' model$beta
tune.fit <- function(x, y, family = c("gaussian", "binomial", "poisson", "cox", "multinom"), penalty = c("SCAD", "MCP", "lasso", "aenet", "msaenet", "enet"), concavity.parameter = switch(penalty, SCAD = 3.7, 3), tune = c("cv", "aic", "bic", "ebic"), nfolds = 10,
                     type.measure = c("deviance", "class", "auc", "mse", "mae"), gamma.ebic = 1, parallel=TRUE, seed=NULL) {
  if(!is.null(seed)){
    set.seed(seed, kind = "L'Ecuyer-CMRG")}
  if (is.null(x) || is.null(y)) {
    stop("The data is missing!")
  }
  
  this.call <- match.call()
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if (class(concavity.parameter) != "numeric") {
    stop("concavity.parameter must be numeric!")
  }
  tune <- match.arg(tune)
  if (class(nfolds) != "numeric") {
    stop("nfolds must be numeric!")
  }
  type.measure <- match.arg(type.measure)
  
  
  if (tune == "cv") {
    
    if (penalty == "lasso") {
      cv.fit <- cv.glmnet(x, y, family = family, type.measure = type.measure, nfolds = nfolds, parallel=parallel)
      coef.beta <- coef(cv.fit, s = "lambda.1se")
      reg.fit <- cv.fit$glmnet.fit
      lambda <- cv.fit$lambda.1se
      lambda.ind <- which(cv.fit$lambda == cv.fit$lambda.1se)
      
    } else if (penalty == "enet"){
      cv.fit = cv.glmnet(x, y, family = family, alpha=0.05, type.measure = type.measure, nfolds = nfolds, parallel=parallel)
      coef.beta = coef(cv.fit, s = "lambda.1se")
      reg.fit = cv.fit$glmnet.fit
      lambda = cv.fit$lambda.1se
      lambda.ind = which(cv.fit$lambda == cv.fit$lambda.1se)
      
    } else if (penalty == "aenet" | penalty == "msaenet") {
      enet <- cv.glmnet(x, y, alpha=0.05, family=family, type.measure = type.measure, nfolds=nfolds, parallel=parallel)
      lambda <- enet$lambda.1se
      coef_init = coef.glmnet(enet, s = lambda)
      if(any(row.names(coef_init)=='(Intercept)')){
        pf=as.vector(pmax(abs(coef_init), .Machine$double.eps)^(-1))[-1]
      } else{
        pf=as.vector(pmax(abs(coef_init), .Machine$double.eps)^(-1))}
      if (penalty == "aenet"){
        if(family == 'cox'){
          cv.fit = Coxnet(x, y, penalty = "Enet", alpha = 0.05, nlambda = 50, nfolds = 10,
                          inzero = FALSE, adaptive = c(TRUE, FALSE), aini=list(wbeta=(pmax(abs(as.numeric(coef_init)), .Machine$double.eps))^(-1)))
          coef.beta = as(as.matrix(cv.fit[['Beta']]), "dgCMatrix")
          rownames(coef.beta) <- colnames(x)
          reg.fit = cv.fit[['fit']]
          lambda = cv.fit[['lambda.max']]
        } else {
          if (family=='gaussian'){
            gcd <- cv.gcdnet(x, y, nfolds=10, lambda2=0.95, standardize=TRUE, method='ls')
            coef_init <- coef(gcd, s = "lambda.1se")
            lambda <- gcd$lambda.1se
            if(any(row.names(coef_init)=='(Intercept)')){
              pf=as.vector(pmax(abs(coef_init), .Machine$double.eps)^(-1))[-1]
            } else{
              pf=as.vector(pmax(abs(coef_init), .Machine$double.eps)^(-1))
            }
            fit_gcdnet <- gcdnet(x, y, lambda=lambda, lambda2=0.95, pf=pf, standardize=TRUE, method='ls')
          } else {
            gcd <- cv.gcdnet(x, y, nfolds=10, lambda2=0.95, standardize=TRUE)
            coef_init <- coef(gcd, s = "lambda.1se")
            lambda = gcd$lambda.1se
            if(any(row.names(coef_init)=='(Intercept)')){
              pf=as.vector(pmax(abs(coef_init), .Machine$double.eps)^(-1))[-1]
            } else{
              pf=as.vector(pmax(abs(coef_init), .Machine$double.eps)^(-1))
            }
            fit_gcdnet <- gcdnet(x, y, lambda=lambda, lambda2=0.95, pf=pf, standardize=TRUE)
          }
          coef.beta = coef(fit_gcdnet, s = "lambda.1se")
          reg.fit = fit_gcdnet$gcdnet.fit
          lambda = fit_gcdnet$lambda.1se
        }
      } else if(penalty == 'msaenet'){
        cv.fit = aenet(x, y, family = family, init = "ridge", alphas = 0.05, penalty.factor.init=pf,
                       rule = "lambda.1se", parallel = parallel, verbose = TRUE)
        coef.beta = cv.fit[['beta']]
        reg.fit = cv.fit[['model']]
        lambda =cv.fit[['best.lambda.aenet']]
        
        # SCAD and MCP penalties    
      }} else if (family != 'cox'){
        cv.fit <- cv.ncvreg(x, y, family = family, penalty = penalty, gamma = concavity.parameter, nfolds = nfolds)
        cv.1se.ind <- min(which(cv.fit$cve < cv.fit$cve[cv.fit$min] + cv.fit$cvse[cv.fit$min]))
        coef.beta <- cv.fit$fit$beta[, cv.1se.ind] # extract coefficients at a single value of lambda, including the intercept
        reg.fit <- cv.fit$fit
        
        lambda <- cv.fit$lambda[cv.1se.ind]
        lambda.ind <- cv.1se.ind
      } else {
        cv.fit <- cv.ncvsurv(x, y, family = family, penalty = penalty, gamma = concavity.parameter, nfolds = nfolds)
        cv.1se.ind <- min(which(cv.fit$cve < cv.fit$cve[cv.fit$min] + cv.fit$cvse[cv.fit$min]))
        coef.beta <- cv.fit$fit$beta[, cv.1se.ind] # extract coefficients at a single value of lambda
        reg.fit <- cv.fit$fit
        
        lambda <- cv.fit$lambda[cv.1se.ind]
        lambda.ind <- cv.1se.ind
      }
  } else {
    n <- nrow(x)
    if (penalty == "lasso") {
      reg.fit <- glmnet(x, y, family = family)
      if(family != 'multinom'){
        coef.beta <- rbind(reg.fit$a0, as.matrix(reg.fit$beta)) # extract coefficients at all values of lambda,  including the intercept
      } else{
        coef.beta = lapply(1:length(reg.fit$beta), f<-function(i){rbind(reg.fit$a0[i,], as.matrix(reg.fit$beta[[i]]))} )
      }
      dev <- deviance(reg.fit)
      reg.df <- reg.fit$df
    } else if (penalty == "enet"){
      reg.fit <- glmnet(x, y, family = family, alpha=0.05)
      if(family != 'multinom'){
        coef.beta <- rbind(reg.fit$a0, as.matrix(reg.fit$beta)) # extract coefficients at all values of lambda,  including the intercept
      } else{
        coef.beta = lapply(1:length(reg.fit$beta), f<-function(i){rbind(reg.fit$a0[i,], as.matrix(reg.fit$beta[[i]]))} )
      }
      dev <- deviance(reg.fit)
      reg.df <- reg.fit$df
    }
    else {
      if (family != "cox") {
        reg.fit <- ncvreg(x, y, family = family, penalty = penalty, gamma = concavity.parameter)
        coef.beta <- reg.fit$beta # extract coefficients at all values of lambda, including the intercept
        dev <- loglik(x, y, coef.beta, family = family)
        reg.df <- getdf(coef.beta[-1, , drop = FALSE])
      } else {
        reg.fit <- ncvsurv(x, y, family = family, penalty = penalty, gamma = concavity.parameter)
        coef.beta <- reg.fit$beta # extract coefficients at all values of lambda, including the intercept
        dev <- 2 * reg.fit$loss
        reg.df <- getdf(coef.beta)
      }
    }
    
    if (tune == "aic") {
      obj <- dev + 2 * reg.df
    }
    if (tune == "bic") {
      obj <- dev + log(n) * reg.df
    }
    if (tune == "ebic") {
      obj <- dev + log(n) * reg.df + 2 * gamma.ebic * log(choose(dim(x)[2], reg.df))
    }
    lambda.ind <- which.min(obj)
    if(family != 'multinom'){
      coef.beta <- coef.beta[, lambda.ind]  
    } else{
      coef.beta = sapply(1:length(coef.beta), f<-function(i){coef.beta[[i]][, lambda.ind]} )
    }
    lambda <- reg.fit$lambda[lambda.ind]
  }
  if (family == 'multinom'){
    a0 <- coef.beta[1,]
    coef.beta <- coef.beta[-1,]  
    ix <- which(apply(abs(coef.beta),1,sum) != 0)
    beta <- coef.beta[ix,]
  } else if (family != "cox") {
    a0 <- coef.beta[1]
    coef.beta <- coef.beta[-1]
    ix <- which(coef.beta != 0)
    beta <- coef.beta[ix]
  } else {
    a0 <- NULL
    coef.beta <- as.vector(coef.beta)
    ix <- which(coef.beta != 0)
    beta <- coef.beta[ix]
  }
  
  return(list(ix = ix, a0 = a0, beta = beta, fit = reg.fit, lambda = lambda))
}
