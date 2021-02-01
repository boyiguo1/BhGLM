get.na.action <- function(na.action) {
  ## get the name of the na.action whether function or text string.
  ## avoids deparse(substitute(na.action)) which is easily broken by 
  ## nested calls.
  if (is.character(na.action)) {
    if (na.action%in%c("na.omit","na.exclude","na.pass","na.fail")) 
      return(na.action) else stop("unrecognised na.action")
  } 
  if (!is.function(na.action)) stop("na.action not character or function")
  a <- try(na.action(c(0,NA)),silent=TRUE)
  if (inherits(a,"try-error")) return("na.fail")
  if (inherits((attr(a,"na.action")),"omit")) return("na.omit")
  if (inherits((attr(a,"na.action")),"exclude")) return("na.exclude")
  return("na.pass")
} ## get.na.action