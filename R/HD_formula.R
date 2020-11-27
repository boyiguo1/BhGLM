# Create spline formula for high-dimension applications
#' Title
#'
#' @param formula an object of class "\link[stats]{formula}": a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’.
#' @param data 
#' @param spl_df 
#' @param rm_overlap 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
#' # spl_df template
#' spline_df <- tribble(
#'  ~Var, ~Func, ~Args,
#'  "x0", "s", "k=1, bs='cr'",
#'  "x1", "s", "k=2, bs='cr'",
#'  "x2", "", "k=2", # ignore Args 
#' )
create_HD_formula <- function(formula, data, spl_df, rm_overlap = TRUE, verbose=T){
  if(missing(formula) & missing(data))
    stop("Either formula or data must be provided")
  
  if(missing(formula)){
    if(!is.data.frame(data)) data <- data.frame(data)
    formula = DF2formula(data)
  }
  
  if(!missing(data)){
    warning("Both formula and dat provided, dat is ignored.")
    warning("Please consider use the function DF2formula and update to improve your formula.")
    # formula[[3]] <- paste(". +", formula[[3]])
    # if(!is.data.frame(dat)) dat <- data.frame(dat)
    # formula <- update(DF2formula(dat),
    #                   formula)  #TODO:update DF2formula given
    
  }
    
  # Manipulate spl_df
  sp_trm <-  spl_df %>% 
    dplyr::filter(Func!="")  %>% # Removing unnecessary terms
    glue::glue_data("{Func}( {Var}{ifelse(is.na(Args)||Args=='', '', paste0(',', Args))})") %>% 
    paste(collapse  = " + ")
  sp_trm <- paste0("~ . + " , sp_trm)
  
  # Adding Spline Terms
  ret <- update(formula, sp_trm)
  
  if(verbose){
    cat("Create formula:\n")
    print(ret)   
  }

  return(ret)
}

