% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmlasso.r
\name{bmlasso}
\alias{bmlasso}
\title{Bayesian Spike-and-Slab Lasso Models}
\usage{
bmlasso(
  x,
  y,
  family = c("gaussian", "binomial", "poisson", "cox"),
  offset = NULL,
  epsilon = 1e-04,
  maxit = 50,
  init = NULL,
  ss = c(0.04, 0.5),
  group = NULL,
  Warning = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{x}{input matrix, of dimension nobs x nvars; each row is an observation vector.}

\item{y}{response variable. Quantitative for family="gaussian", or family="poisson" (non-negative counts). For family="gaussian", y is always been standardized. For family="binomial", y should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). For family="cox", y should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' indicating right censored. The function Surv() in package \bold{survival} produces such a matrix.}

\item{family}{}

\item{offset}{}

\item{epsilon}{}

\item{maxit}{}

\item{init}{}

\item{ss}{}

\item{group}{}

\item{Warning}{}

\item{verbose}{}
}
\value{
This function returns all outputs from the function \code{\link{glmnet}}, and some other values used in Bayesian hierarchical models.
}
\description{
This function is to set up Bayesian 
GLMs or Cox survival models with spike-and-slab 
mixture double-exponential prior (called the Bayesian
spike-and-slab mixture lasso), and to fit the model
 by incorporating EM steps into the fast coordinate
  descent algorithm.
}
\examples{
library(BhGLM)
library(survival)
library(glmnet)


N = 1000
K = 100
x = sim.x(n=N, m=K, corr=0.6) # simulate correlated continuous variables  
h = rep(0.1, 4) # assign four non-zero main effects to have the assumed heritabilty 
nz = as.integer(seq(5, K, by=K/length(h))); nz
yy = sim.y(x=x[, nz], mu=0, herit=h, p.neg=0.5, sigma=1.6) # simulate responses
yy$coefs

# y = yy$y.normal; fam = "gaussian"; y = scale(y)
# y = yy$y.ordinal; fam = "binomial"
y = yy$y.surv; fam = "cox" 

group = NULL
#group = rep(0, 21)
#for(j in 1:length(group)) group[j] = (j-1) * K/(length(group)-1)

# lasso and mixture lasso

f1 = glmNet(x, y, family = fam, ncv = 1) 

ps = f1$prior.scale; ps 
ss = c(ps, 0.5)
f2 = bmlasso(x, y, family = fam, ss = ss, group = group)

par(mfrow = c(1, 2), mar = c(3, 4, 4, 4))
gap = 10
plot.bh(coefs = f1$coef, threshold = f1$df, gap = gap, main = "lasso") 
plot.bh(coefs = f2$coef, threshold = f2$df, gap = gap, main = "mixture lasso") 

}
\references{
{
Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for Generalized Linear Models via Coordinate Descent. J Stat Softw 33, 1-22.

Simon, N., Friedman, J., Hastie, T. & Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. Journal of Statistical Software 39, 1-13.

  Zaixiang Tang, Yueping Shen, Xinyan Zhang, Nengjun Yi (2017) The Spike-and-Slab Lasso Generalized Linear Models for Prediction and Associated Genes Detection. Genetics 205, 77-88.
  
  Zaixiang Tang, Yueping Shen, Xinyan Zhang, Nengjun Yi (2017) The Spike-and-Slab Lasso Cox Models for Survival Prediction and Associated Genes Detection. Bioinformatics, 33(18), 2799-2807.

  Zaixiang Tang, Yueping Shen, Yan Li, Xinyan Zhang, Jia Wen, Chen'ao Qian, Wenzhuo Zhuang, Xinghua Shi, and Nengjun Yi (2018) Group Spike-and-Slab Lasso Generalized Linear Models for Disease Prediction and Associated Genes Detection by Incorporating Pathway Information. Bioinformatics 34(6): 901-910.

Zaixiang Tang, Yueping Shen, Shu-Feng Lei, Xinyan Zhang, Zixuan Yi, Boyi Guo, Jake Chen, and Nengjun Yi (2019) Gsslasso Cox: a fast and efficient pathway-based framework for predicting survival and detecting associated genes. BMC Bioinformatics 20(94). 
}
}
\seealso{
{
  \code{\link{glmnet}}, \code{\link{glmNet}}, \code{\link{bglm}}, \code{\link{bcoxph}}
}
}
\author{
{
 Nengjun Yi, nyi@uab.edu
}
}
