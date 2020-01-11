

#' Bayesian Spike-and-Slab Lasso Models
#'
#' This function is to set up Bayesian 
#'GLMs or Cox survival models with spike-and-slab 
#'mixture double-exponential prior (called the Bayesian
#' spike-and-slab mixture lasso), and to fit the model
#'  by incorporating EM steps into the fast coordinate
#'   descent algorithm. 
#'
#'
#' @param x input matrix, of dimension nobs x nvars; each row is an observation vector.
#' @param y response variable. Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="gaussian", y is always been standardized. For family="binomial", y should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For family="cox", y should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' indicating right censored. The function Surv() in package \bold{survival} produces such a matrix. 
#' @param family 
#' @param offset 
#' @param epsilon 
#' @param maxit 
#' @param init 
#' @param ss 
#' @param group 
#' @param Warning 
#' @param verbose 
#'
#' @return This function returns all outputs from the function \code{\link{glmnet}}, and some other values used in Bayesian hierarchical models.
#' @export
#' @import glmnet
#' @importFrom stats sd
#' @importFrom survival is.Surv
#'
#' @references {
#' Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for Generalized Linear Models via Coordinate Descent. J Stat Softw 33, 1-22.
#'
#' Simon, N., Friedman, J., Hastie, T. & Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. Journal of Statistical Software 39, 1-13.
#' 
#'   Zaixiang Tang, Yueping Shen, Xinyan Zhang, Nengjun Yi (2017) The Spike-and-Slab Lasso Generalized Linear Models for Prediction and Associated Genes Detection. Genetics 205, 77-88.
#'   
#'   Zaixiang Tang, Yueping Shen, Xinyan Zhang, Nengjun Yi (2017) The Spike-and-Slab Lasso Cox Models for Survival Prediction and Associated Genes Detection. Bioinformatics, 33(18), 2799-2807.
#' 
#'   Zaixiang Tang, Yueping Shen, Yan Li, Xinyan Zhang, Jia Wen, Chen'ao Qian, Wenzhuo Zhuang, Xinghua Shi, and Nengjun Yi (2018) Group Spike-and-Slab Lasso Generalized Linear Models for Disease Prediction and Associated Genes Detection by Incorporating Pathway Information. Bioinformatics 34(6): 901-910.
#' 
#' Zaixiang Tang, Yueping Shen, Shu-Feng Lei, Xinyan Zhang, Zixuan Yi, Boyi Guo, Jake Chen, and Nengjun Yi (2019) Gsslasso Cox: a fast and efficient pathway-based framework for predicting survival and detecting associated genes. BMC Bioinformatics 20(94). 
#' }
#' 
#' @author{
#'  Nengjun Yi, nyi@uab.edu
#' }
#' 
#' @seealso{
#'   \code{\link{glmnet}}, \code{\link{glmNet}}, \code{\link{bglm}}, \code{\link{bcoxph}}
#' }
#'
#' @examples
#' library(BhGLM)
#' library(survival)
#' library(glmnet)
#' 
#' 
#' N = 1000
#' K = 100
#' x = sim.x(n=N, m=K, corr=0.6) # simulate correlated continuous variables  
#' h = rep(0.1, 4) # assign four non-zero main effects to have the assumed heritabilty 
#' nz = as.integer(seq(5, K, by=K/length(h))); nz
#' yy = sim.y(x=x[, nz], mu=0, herit=h, p.neg=0.5, sigma=1.6) # simulate responses
#' yy$coefs
#' 
#' # y = yy$y.normal; fam = "gaussian"; y = scale(y)
#' # y = yy$y.ordinal; fam = "binomial"
#' y = yy$y.surv; fam = "cox" 
#' 
#' group = NULL
#' #group = rep(0, 21)
#' #for(j in 1:length(group)) group[j] = (j-1) * K/(length(group)-1)
#' 
#' # lasso and mixture lasso
#' 
#' f1 = glmNet(x, y, family = fam, ncv = 1) 
#' 
#' ps = f1$prior.scale; ps 
#' ss = c(ps, 0.5)
#' f2 = bmlasso(x, y, family = fam, ss = ss, group = group)
#' 
#' par(mfrow = c(1, 2), mar = c(3, 4, 4, 4))
#' gap = 10
#' plot.bh(coefs = f1$coef, threshold = f1$df, gap = gap, main = "lasso") 
#' plot.bh(coefs = f2$coef, threshold = f2$df, gap = gap, main = "mixture lasso") 
#' 
bmlasso <- function(x, y, family = c("gaussian", "binomial", "poisson", "cox"), offset = NULL,
                    epsilon = 1e-04, maxit = 50, init = NULL, 
                    ss = c(0.04, 0.5), group = NULL, Warning = FALSE, verbose = FALSE) 
{
  # if (!requireNamespace("glmnet")) install.packages("glmnet")
  # require(glmnet)
  start.time <- Sys.time()
  call <- match.call()
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste("x", 1:ncol(x), sep = "")
  nobs <- nrow(x)
  if (NROW(y) != nobs) stop("nobs of 'x' and 'y' are different")
  inc <- apply(cbind(y, x), 1, function(z) !any(is.na(z)))
  if (!is.null(offset)) {
    if (length(offset) != nobs) stop("nobs of 'x' and 'offset' are different")
    inc <- apply(cbind(y, x, offset), 1, function(z) !any(is.na(z)))
  }
  y <- y[inc]
  x <- x[inc,]
  offset <- offset[inc]
  family <- family[1]
  if (family == "cox")  
    if (!is.Surv(y)) stop("'y' should be a 'Surv' object")
  if (family == "gaussian") y <- (y - mean(y))/sd(y)
  if (!is.null(init) & length(init) != ncol(x)) stop("give an initial value to each coefficient (not intercept)")

  f <- bmlasso.fit(x = x, y = y, family = family, offset = offset, epsilon = epsilon, maxit = maxit, init = init,
                   group = group, ss = ss, Warning = Warning)
  
  f$call <- call
  if (family == "cox") class(f) <- c(class(f), "bmlasso", "COXPH") 
  else class(f) <- c(class(f), "bmlasso", "GLM")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"), 3)
  if (verbose){
    cat("EM Coordinate Decent Iterations:", f$iter, "\n")
    cat("Computational time:", minutes, "minutes \n")
  }

  return(f)
}

# ******************************************************************************

bmlasso.fit <- function(x, y, family = "gaussian", offset = NULL, epsilon = 1e-04, maxit = 50, 
                        init = rep(0, ncol(x)), group = NULL, ss = c(0.04, 0.5), 
                        Warning = FALSE)
{  
  ss <- sort(ss)
  ss <- ifelse(ss <= 0, 0.001, ss)
  prior.scale <- ss[length(ss)]  # used for ungrouped coefficients 

  if (family == "cox") intercept <- FALSE
  else intercept <- TRUE
  x0 <- x
  if (intercept) x0 <- cbind(1, x)
  d <- prepare(x = x0, intercept = intercept, prior.mean = 0, prior.sd = 1, prior.scale = prior.scale, 
               prior.df = 1, group = group)
  x <- d$x
  prior.scale <- d$prior.scale 
  group <- d$group
  group.vars <- d$group.vars
  ungroup.vars <- d$ungroup.vars
  prior.scale <- prior.scale / autoscale(x, min.x.sd=1e-04)
  if (intercept){
    x <- x[, -1]
    prior.scale <- prior.scale[-1]
  }
  
  if (length(ss) != 2) stop("ss should have two positive values")
  gvars <- unlist(group.vars)
  theta <- p <- rep(0.5, length(gvars))
  names(theta) <- names(p) <- gvars 
  
  if (is.null(init)) {
    for (k in 1:5) {
      ps <- ss[1] + (k - 1) * 0.01
      if (family == "cox") ps <- min(ss[1] + (k - 1) * 0.01, 0.08)
      f <- glmnet(x=x, y=y, family=family, offset=offset, alpha=0.95, 
                  lambda=1/(nrow(x) * ps), standardize=TRUE)
      b <- as.numeric(f$beta)
      if (any(b != 0)) break
    }
  }
  else b <- as.numeric(init)
  names(b) <- colnames(x)
  b <- ifelse(b == 0, 0.001, b)
  init <- b
  
  devold <- 0
  conv <- FALSE
  for (iter in 1:maxit){
    
    out <- update.scale.p(b0=b[gvars], ss=ss, theta=theta)
    prior.scale[gvars] <- out[[1]]   
    p <- out[[2]]
    if (!is.matrix(group))
      theta <- update.ptheta.group(group.vars=group.vars, p=p)
    else theta <- update.ptheta.network(theta=theta, p=p, w=group)
    
    Pf <- 1/(prior.scale + 1e-10)
    f <- glmnet(x = x, y = y, family = family, offset = offset, alpha = 1, 
                penalty.factor = Pf, lambda = sum(Pf)/(nrow(x) * ncol(x)), standardize = FALSE)
    
    b <- as.numeric(f$beta) #/sqrt(dispersion)
    names(b) <- colnames(x)
    dev <- deviance(f)
  
    if(abs(dev - devold)/(0.1 + abs(dev)) < epsilon & iter > 5) {
      conv <- TRUE
      break
    }
    else devold <- dev
  }
  if (Warning & !conv) warning("algorithm did not converge", call. = FALSE)
  
  f$x <- x
  f$y <- y
  f$family <- family
  f$ss <- ss
  
  f$coefficients <- as.numeric(coef(f))
  names(f$coefficients) <- rownames(coef(f))
  f$linear.predictors <- predict(f, newx = x, type = "link", offset = offset)
  if (family == "gaussian")
    f$dispersion <- bglm(y ~ f$linear.predictors-1, start=1, prior=De(1,0), verbose=FALSE)$dispersion 
  
  f$iter <- iter
  f$prior.scale <- prior.scale
  f$penalty.factor <- Pf
  f$group <- group
  f$group.vars <- group.vars 
  f$ungroup.vars <- ungroup.vars
  f$p <- p
  f$ptheta <- theta
  f$init <- init
  f$aic <- deviance(f) + 2 * f$df
  f$offset <- offset
  
  return(f)
}

#*******************************************************************************

