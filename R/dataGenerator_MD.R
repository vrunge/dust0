

####################################################
#############    dataGenerator_MD   ################
####################################################

#' dataGenerator_MD
#'
#' @description Generating copies of univariate time series for multiple change-point detection based on uni-parametric models of the exponential family
#' @param chpts a vector of increasing change-point indices (the last value is data length)
#' @param parameters dataframe of successive segment parameters (each column = parameters for one time-series) (as many rows as values in \code{chpts} vector)
#' @param sdNoise (type \code{"gauss"}) standard deviation for the noise parameter. If a vector, value at position i for i-th time series.
#' @param gamma (type \code{"gauss"}) vector or dataframe (each column = parameters for one time-series) of numbers between 0 and 1 : the coefficient of the exponential decay. By default = 1 for piecewise constant signals. If one value, it is used for all segments. Otherwise we need as many values as in \code{chpts} vector.
#' @param nbTrials  (type \code{"binom"}) number of trials. If a vector, value at position i for i-th time series.
#' @param nbSuccess (type \code{"negbin"}) number of successes. If a vector, value at position i for i-th time series.
#' @param type the model: \code{"gauss"}, \code{"exp"}, \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}
#' @return a multivariate time series (copies of univariate time-series with the same model type) following the chosen model type and parameters
#' @examples
#' dataGenerator_MD(chpts = c(50,100),
#'                  parameters = data.frame(ts1 = c(10,10), ts2 = c(5,20)),
#'                  gamma = data.frame(ts1 = c(0.9,0.8), ts2 = c(0.8,0.9)),
#'                  sdNoise = 0.2, type = "gauss")
#' dataGenerator_MD(chpts = c(50,100,150),
#'                  parameters = data.frame(ts1 = c(4,0,2), ts2 = c(5,2,-2),ts3 = c(0,1,0)),
#'                  sdNoise = 0.4, type = "gauss")
#' dataGenerator_MD(chpts = c(50,100,150),
#'                  parameters = data.frame(ts1 = c(4,10,2), ts2 = c(5,2,2),ts3 = c(10,4,3)),
#'                  type = "poisson")
#' dataGenerator_MD(chpts = c(50,100,150),
#'      parameters = data.frame(ts1 = c(0.4,0.3,0.5), ts2 = c(0.5,0.6,0.2),ts3 = c(0.1,0.4,0.6)),
#'      nbTrials = c(30,10,100), type = "binom")
dataGenerator_MD <- function(chpts = 100,
                             parameters = data.frame(ts1 = 0, ts2 = 0),
                             sdNoise = 1,
                             gamma = 1,
                             nbTrials = 10,
                             nbSuccess = 10,
                             type = "gauss")
{
  p <- ncol(parameters)
  if(p == 1){stop('use the function dataGenerator_1D')}

  ##################################################
  ### replicate parameters for multivariate case ###
  ##################################################

  if(length(sdNoise) == 1){sdNoise <- rep(sdNoise, p)}
  if(!is.data.frame(gamma)){gamma <- data.frame(matrix(gamma, nrow = length(gamma), ncol = p))}
  if(length(nbTrials) == 1){nbTrials <- rep(nbTrials, p)}
  if(length(nbSuccess) == 1){nbSuccess <- rep(nbSuccess, p)}

  ############  data generation   ############

  res <- matrix(NA, nrow = p, ncol = chpts[length(chpts)])
  for(i in 1:p)
  {
    res[i,] <- dataGenerator_1D(chpts = chpts,
                                parameters = parameters[,i],
                                sdNoise = sdNoise[i],
                                gamma = gamma[,i],
                                nbTrials = nbTrials[i],
                                nbSuccess = nbSuccess[i],
                                type = type)
  }
  return(res)
}
