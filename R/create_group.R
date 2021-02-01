#' Creating a numeric vector to fit bglm_spline
#'
#' @param G The mgcv pre-fit object that contains smooth spline information
#' @param drop.intercept 
#'
#' @return
#' @export
#'
#' @examples
#' N <- 1000
#' p <- 10
#' 
#' set.seed(1)
#' 
#' dat <- data.frame(sim_Bai_logistic(N, p)$dat) 
#' G <- gam(y ~ s(x1, bs = "ps", k=13)+s(x2, bs = "ps", k=15)+s(x3, bs = "ps", k=15)+s(x4, bs = "ps", k=15)+
#' s(x5, bs = "ps", k=15)+s(x6, bs = "ps", k=15)+s(x7, bs = "ps", k=15)+s(x8, bs = "ps", k=15)+
#'  s(x9, bs = "ps", k=15)+s(x10, bs = "ps", k=15), 
#' data=dat %>% data.frame(), family=binomial,drop.intercept = TRUE, fit = FALSE)
#' 
#' colnames(G$X) <- G$term.names
#' bgam_mdl <- bglm_spline(G$y~G$X-1, family = G$family, prior=mt(mean=0, s0=0.04, s1=1, df=Inf),
#'                         group = create_group(G))

create_group <- function(G, drop.intercept=FALSE)
{
  m <- length(G$smooth)
  
  group <- c(G$smooth[[1]]$first.para-1,
             map_dbl(G$smooth, .f = function(x){x$last.para}))
  
  
  if(drop.intercept == TRUE)
    group <- group - 1
  # TODO: consider if drop.intercept is needed. in this case, drop.intercept should be aligned with if the new formula have intercept, which possibly alter the position of the groups
  return(group)
}