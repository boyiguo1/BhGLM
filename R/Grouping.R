
# What does Grouping do
# Check if group is a matrix
# # if yes, calculate new group weights for future use by |x_{ij}|/(\sum\limits^n_{t=1}|x_{it}|-|x_{ii}|)
# When group is not matrix
# Check if group is a list,
# If not, doing the old fashion thing where
# 



Grouping <- function(all.var, group) 
{ 
  n.vars <- length(all.var)
  group.vars <- list()
  
  if (is.matrix(group))
  {
    if (nrow(group)!=ncol(group) | ncol(group)>n.vars) 
      stop("wrong dimension for 'group'")
    if (any(rownames(group)!=colnames(group)))
      stop("rownames should be the same as colnames")
    if (any(!colnames(group)%in%all.var))
      stop("variabe names in 'group' not in the model predictors")
    group.vars <- colnames(group)
    group <- abs(group)
    wcol <- rowSums(group) - diag(group)
    group <- group/wcol
  }
  else{
    if (is.list(group)) group.vars <- group
    else
    {
      if (is.numeric(group) & length(group)>1) { 
        group <- sort(group)  
        if (group[length(group)] > n.vars) stop("wrong grouping")
      }
      if (is.numeric(group) & length(group)==1)
        group <- as.integer(seq(0, n.vars, length.out = n.vars/group + 1))
      if (is.null(group)) group <- c(0, n.vars)
      group <- unique(group)
      for (j in 1:(length(group) - 1))
        group.vars[[j]] <- all.var[(group[j] + 1):group[j + 1]]
    }
  }
  all.group.vars <- unique(unlist(group.vars))
  
  if (length(all.group.vars) == n.vars) ungroup.vars <- NULL
  else ungroup.vars <- all.var[which(!all.var %in% all.group.vars)]
  
  group.new <- c(length(ungroup.vars), length(ungroup.vars) + cumsum(lapply(group.vars, length)))
  var.new <- c(ungroup.vars, unlist(group.vars))
  
  list(group=group, # the new "distance" matrix used in the network analysis, other wise, old thing
       group.vars=group.vars, # viables that appear in the group
       ungroup.vars=ungroup.vars, # variable that didn't appear in the group parameter
       group.new=group.new, # no idea what does this do
       var.new=var.new   # all variables ordered by ungrouped variables, and grouped variables
  ) 
}