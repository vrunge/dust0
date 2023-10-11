
#' OP_R
#'
#' @description OP algorithm for univariate time-series (with different possible data models)
#' @param data a vector of data
#' @param penalty penalty value (non-negative)
#' @param type type of cost to use: gauss, poisson, exp
#' @return a list with the change-point elements (each last index of each segment) and a vector `nb` saving the number of non-pruned elements at each iteration
#' @examples
#' OP_R(dataGenerator_1D(chpts = c(50,200,400), parameters = c(0,1,0), type = "gauss"), 2*log(400))
OP_R <- function(data, penalty, type = "gauss")
{
  ##########  ##########  ##########  ##########  ##########

  if(!is.vector(data)){stop('data is not a vector')}
  if(length(data) <= 1){stop('no data to segment')}
  if(penalty < 0){stop('penalty must be non negative')}
  allowed.types <- c("gauss", "poisson", "exp")
  if(!type %in% allowed.types){stop('type must be one of the list: ', paste(allowed.types, collapse=", "))}

  ##########  ##########  ##########  ##########  ##########

  if(type == "gauss"){res <- OP_R_gauss(data, penalty)}
  #if(type == "poisson"){res <- dust_R_poisson(data, penalty)}
  #if(type == "exp"){res <- dust_R_exp(data, penalty)}

  return(res)
}


##################################################
##################################################
##################################################



OP_R_gauss <- function(data, penalty)
{
  #########
  ###
  ### DATA preprocessing
  ###
  n <- length(data)
  cumData <- cumsum(c(0, data))
  cumData2 <- cumsum(c(0, data^2))

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
      eval <- costQ[shift(k-1)] + (t-k+1)*eval_var(cumData, cumData2, k, t) + penalty
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t)] <- min_temp
    cp[shift(t)] <- index - 1

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
  return(list(changepoints = changepoints[-1], costQ = costQ))
}






