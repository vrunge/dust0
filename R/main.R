

#' dust_R
#'
#' @description DUST algorithm for univariate time-series (with different possible data models)
#' @param data a vector of data
#' @param penalty penalty value (non-negative)
#' @param type type of cost to use: gauss, poisson, exp
#' @return a list with the change-point elements (each last index of each segment) and a vector `nb` saving the number of non-pruned elements at each iteration
#' @examples
#' dust_R(dataGenerator_1D(chpts = c(50,200,400), parameters = c(0,1,0), type = "gauss"), 2*log(400))
dust_R <- function(data, penalty, type = "gauss")
{
  ##########  ##########  ##########  ##########  ##########

  if(!is.vector(data)){stop('data is not a vector')}
  if(length(data) <= 1){stop('no data to segment')}
  if(penalty < 0){stop('penalty must be non negative')}
  allowed.types <- c("gauss", "poisson", "exp")
  if(!type %in% allowed.types){stop('type must be one of the list: ', paste(allowed.types, collapse=", "))}

  ##########  ##########  ##########  ##########  ##########

  if(type == "gauss"){res <- dust_R_gauss(data, penalty)}
  if(type == "poisson"){res <- dust_R_poisson(data, penalty)}
  if(type == "exp"){res <- dust_R_exp(data, penalty)}

  return(res)
}


##################################################
##################################################
##################################################





