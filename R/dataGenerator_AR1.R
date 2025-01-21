

###############################################################################################
## NON CLOSED COST

##################################################
#############  dataGenerator_AR1  ################
##################################################

#' dataGenerator_AR1
#'
#' @description Generating time series for multiple change-point detection with AR1 model
#' @param chpts a vector of increasing change-point indices (the last value is data length)
#' @param means vector of successive segment parameters (as many parameters as values in \code{chpts} vector)
#' @param sdNoise standard deviation for the noise parameter (type \code{"gauss"})
#' @param rho vector of numbers between 0 and 1 : the coefficient of the exponential decay (type \code{"gauss"}). By default = 1 for piecewise constant signals. If one value, it is used for all segments. Otherwise we need as many values as in \code{chpts} vector.
#' @return a time series following the AR(1) model type
#' @examples
#' dataGenerator_AR1(chpts = c(50,100), means = c(0,1), sdNoise = 0.2, rho = c(0.9,0.7))
dataGenerator_AR1 <- function(chpts = 100,
                              means = 0,
                              sdNoise = 1,
                              rho = 0.7)
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(chpts)){stop('chpts values are not all numeric')}
  if(!all(chpts > 0)){stop('chpts values are not all positives')}
  if(is.unsorted(chpts, strictly = TRUE)){stop('chpts should be a strictly increasing vector of change-point positions (indices)')}

  if(!is.numeric(means)){stop('means values are not all numeric')}
  if(length(chpts) != length(means)){stop('chpts and means vectors are of different size')}

  if(!is.numeric(rho)){stop('rho values are not all numeric')}
  {if(any(rho > 1) || any(rho < 0)){stop('The vector rho should contain values between 0 and 1')}}
  if(length(chpts) != length(rho)){stop('chpts and rho vectors are of different size')}

  #############  #############  #############
  ############ data generation   ############
  #############  #############  #############

  n <- chpts[length(chpts)]
  repetition <- c(chpts[1], diff(chpts))
  mu <- rep(means, repetition)


  y <- rnorm(n, mean = mu, sd = sdNoise)

  return(y)
}



