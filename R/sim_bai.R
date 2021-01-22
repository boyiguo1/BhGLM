sim_Bai_logistic <- function(n, p){
  if(p < 4){
    p <- 4
    warning("There must be more than 4 variables. Changing p to 4")
  }
  
  x <- MASS::mvrnorm(n=n, mu= rep(0, p), Sigma = diag(p))
  lp <- 5*sin(2*pi*x[,1]) - 4*cos(2*pi*x[,2]-0.5) + 6*(x[,3]-0.5)-5*(x[,4]^2-0.3)
  theta <- exp(lp)/(1+exp(lp))
  y <- rbinom(n, 1, theta)
  
  colnames(x) <- paste0("x", 1:p)
  
  return(list(dat=cbind(x,y), lp=lp, theta = theta))
}
