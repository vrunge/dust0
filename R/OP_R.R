
#' OP_R
#'
#' @description OP algorithm for univariate time-series (with different possible data models)
#' @param data a vector of data (univariate)
#' @param penalty penalty value (non-negative)
#' @param type type of cost to use: \code{"gauss"}, \code{"exp"}, \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}
#' @return a list with the change-point elements (each last index of each segment) and a vector \code{costQ} saving the optimal cost at each time step
#' @examples
#'  OP_R(dataGenerator_1D(chpts = c(50,100), c(0,1), sdNoise = 1, type = "gauss"), log(100))
#'  OP_R(dataGenerator_1D(chpts = c(50,100), c(1,7), type = "exp"), log(100))
#'  OP_R(dataGenerator_1D(chpts = c(50,100), c(3,10), type = "poisson"), 10*log(100))
#'  OP_R(dataGenerator_1D(chpts = c(50,100), c(0.7,0.3), type = "geom"), 5*log(100))
#'  OP_R(dataGenerator_1D(chpts = c(50,100), c(0.7,0.2), type = "bern"), log(100))
#'  OP_R(dataGenerator_1D(chpts = c(50,100), c(0.7, 0.3), nbTrials = 5, type = "binom"), 5*log(100))
#'  OP_R(dataGenerator_1D(chpts = c(50,100), c(0.4,0.7), nbSuccess = 10, type = "negbin"), 50*log(100))
OP_R <- function(data, penalty, type = "gauss")
{
  ##########  ##########  ##########  ##########  ##########

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


  #########
  ###
  ### update rule Dynamic Programming
  ###
  for(t in 1:n) # at t, transform Q_{t-1} into Q_{t}
  {
    indexSet <- c(indexSet, t) #add new test point

    min_temp <- Inf

    for(k in indexSet)
    {
      eval <- min_cost(A, B, S, k, t+1, costQ[k] + penalty)
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[t+1] <- min_temp
    cp[t+1] <- index-1
  }

  #########
  ###
  ### backtracking step
  ###
  changepoints <- n # vector of change-point to build
  current <- n

  while(changepoints[1] > 0)
  {
    pointval <- cp[shift(current)] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  return(list(changepoints = changepoints[-1], costQ = costQ[-1]))
}




