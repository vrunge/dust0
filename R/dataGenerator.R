

#' dataGenerator
#' @description Generating data for change-point detection for univariate time series
#' @param chpts a vector of increasing change-point indices
#' @param parameters vector of successive segment parameters
#' @param sdNoise standard deviation for the noise parameter (if available)
#' @param type the model: gauss or poisson
#' @return a vector of simulated time series
#' @examples
#' myData <- dataGenerator(chpts = c(30,100,120), parameter = c(0,5,0))
dataGenerator <- function(chpts = 100,
                          parameters = 0,
                          sdNoise = 1,
                          type = "gauss")
{
  ############
  ### STOP ###
  ############

  if(!is.numeric(chpts)){stop('chpts values are not all numeric')}
  if(is.unsorted(chpts)){stop('chpts should be an increasing vector of change-point positions (indices)')}
  if(length(unique(chpts)) < length(chpts)){stop('chpts is not a strictly increasing sequence')}
  if(!is.numeric(parameters)){stop('parameters values are not all numeric')}
  if(length(chpts) != length(parameters)){stop('chpts and parameters vectors are of different size')}
  if(sdNoise < 0){stop('sdNoise cannot be negative')}
  if(type == "poisson"){if(min(parameters) <= 0){stop('no negative mean allowed for Poisson model')}}

  ############
  ############
  ############

  n <- chpts[length(chpts)]
  repetition <- c(chpts[1], diff(chpts))
  mu <- rep(parameters, repetition)

  if(type == "gauss")
  {
    y <- rnorm(n, mean = mu, sd = sdNoise)
  }
  if(type == "poisson")
  {
    y <- rpois(n = n, lambda = mu)
  }
  return(y)
}







