

#' dust_R
#'
#' @description DUST algorithm for univariate time-series (with different possible data models)
#' @param data a vector of data (univariate)
#' @param penalty penalty value (non-negative)
#' @param type type of cost to use: \code{"gauss"}, \code{"exp"}, \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}
#' @return a list with the change-point elements (each last index of each segment) and a vector `\code{nb} saving the number of non-pruned elements at each iteration
#' @examples
#' dust_R(dataGenerator_1D(chpts = c(50,200,400), parameters = c(0,1,0), type = "gauss"), 2*log(400))
dust_R <- function(data, penalty, type = "gauss")
{

  if(!is.vector(data)){stop('data is not a vector')}
  if(length(data) <= 1){stop('no data to segment')}
  if(penalty < 0){stop('penalty must be non negative')}

  allowed.types <- c("gauss", "exp", "poisson", "geom", "bern", "binom", "negbin")
  if(!type %in% allowed.types){stop('type must be one of: ', paste(allowed.types, collapse=", "))}

  ##########  ##########  ##########  ##########  ##########

  ### preprocessing
  n <- length(data)
  stat <- statistic(type = type)
  S <- c(0, cumsum(stat(data)))

  # loading the type specific functions
  A <- A(type = type)
  B <- B(type = type)

  #########
  ###
  ### INITIALIZATION
  ###
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -penalty

  indexSet <- NULL


  return(NULL)
}






