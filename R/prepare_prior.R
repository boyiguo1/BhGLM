# the function is to prepare the piors to the coefficients
# Prior for intercept: mean -> 0, scale <- same as others, df <- same as others
# Prior for variables: specified in the arguments
# When prior.mean is smaller than the number of variables after accounting for the intercept
# # making up the missing ones with the last one in the prior mean
# When the sample data standard deviation is extremely small, setting very small prior sd
# Setting the prior sd for the intercept very large (1e+10)
# Setting the ungrouped variables' prior mean as 0

prepare <- function(x, intercept, prior.mean, prior.sd, prior.scale, prior.df, group)
{
  x0 <- x
  if (intercept) x0 <- x[, -1, drop = FALSE] 
  g <- Grouping(all.var = colnames(x0), group = group)
  group <- g$group # update the distance matrix 
  group.vars <- g$group.vars # get the group variables, in my case, not necessary
  ungroup.vars <- g$ungroup.vars # not necessary
  covars <- g$ungroup.vars  # not necessary
  
  if (is.list(group)) { # for overlap groups
    if (length(unlist(group)) > length(unique(unlist(group)))) {
      x1 <- as.data.frame(x0)
      x1 <- x1[, c(covars, unlist(group))]
      g <- c(length(ungroup.vars), length(ungroup.vars) + cumsum(lapply(group, length)))
      for (j in 1:(length(group)-1))
        group.vars[[j]] <- colnames(x1[, (g[j]+1):g[j+1]])
      x1 <- as.matrix(x1)
      x <- x1 
      if (intercept) {
        x <- cbind(1, x)
        colnames(x)[1] <- "(Intercept)"
      }
    }
  }
  
  J <- NCOL(x)
  
  if (intercept & J > 1) {
    prior.mean <- c(0, prior.mean)
    prior.scale <- c(prior.scale[1], prior.scale)
    prior.df <- c(prior.df[1], prior.df)
  }
  
  if (length(prior.mean) < J) 
    prior.mean <- c(prior.mean, rep(prior.mean[length(prior.mean)], J - length(prior.mean)) )
  if (length(prior.scale) < J) 
    prior.scale <- c(prior.scale, rep(prior.scale[length(prior.scale)], J - length(prior.scale)) )
  if (length(prior.df) < J) 
    prior.df <- c(prior.df, rep(prior.df[length(prior.df)], J - length(prior.df)) )
  prior.mean <- prior.mean[1:J]
  prior.scale <- prior.scale[1:J]
  prior.df <- prior.df[1:J]
  prior.df <- ifelse(prior.df==Inf, 1e+10, prior.df)
  
  if (is.null(prior.sd)) prior.sd <- prior.scale + 0.2   ## + 0.2 to avoid prior.scale=0
  if (length(prior.sd) < J)  
    prior.sd <- c(prior.sd, rep(prior.sd[length(prior.sd)], J - length(prior.sd)) )
  prior.sd <- prior.sd[1:J]
  sd.x <- apply(x, 2, sd, na.rm=TRUE)
  min.x.sd <- 1e-04
  prior.sd <- ifelse(sd.x < min.x.sd, 1e-04, prior.sd)
  if (intercept) prior.sd[1] <- 1e+10 
  
  names(prior.mean) <- names(prior.scale) <- names(prior.df) <- names(prior.sd) <- colnames(x)
  
  if (intercept) covars <- c(colnames(x)[1], covars)
  if (!is.null(covars)) prior.mean[covars] <- 0
  
  list(x=x, prior.mean=prior.mean, prior.sd=prior.sd, prior.scale=prior.scale, prior.df=prior.df, 
       sd.x=sd.x, min.x.sd=min.x.sd,
       group=group, group.vars=group.vars, ungroup.vars=ungroup.vars)
}  