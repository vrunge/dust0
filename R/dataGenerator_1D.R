
####################################################
#############    dataGenerator_1D   ################
####################################################

#' dataGenerator_1D
#'
#' @description Generating univariate time series for multiple change-point detection based on uni-parametric models of the exponential family
#' @param chpts a vector of increasing change-point indices (the last value is data length)
#' @param parameters vector of successive segment parameters (as many parameters as values in \code{chpts} vector)
#' @param sdNoise  (type \code{"gauss"}) standard deviation for the noise parameter
#' @param gamma (type \code{"gauss"}) vector of numbers between 0 and 1 : the coefficient of the exponential decay. By default = 1 for piecewise constant signals. If one value, it is used for all segments. Otherwise we need as many values as in \code{chpts} vector.
#' @param nbTrials (type \code{"binom"}) number of trials
#' @param nbSuccess (type \code{"negbin"}) number of successes
#' @param type the model: \code{"gauss"}, \code{"exp"}, \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}, \code{"variance"}
#' @return a univariate time series following the chosen model type and parameters
#' @examples
#' dataGenerator_1D(chpts = c(50,100), parameters = c(0,1), sdNoise = 0.2, type = "gauss")
#' dataGenerator_1D(chpts = c(50,100), parameters = c(10,20), gamma = c(0.9,0.95), type = "gauss")
#' dataGenerator_1D(chpts = c(50,100), parameters = c(2,7), type = "exp")
#' dataGenerator_1D(chpts = c(50,100), parameters = c(3,5), type = "poisson")
#' dataGenerator_1D(chpts = c(50,100), parameters = c(0.6,0.3), type = "geom")
#' dataGenerator_1D(chpts = c(50,100), parameters = c(0.7,0.2), type = "bern")
#' dataGenerator_1D(chpts = c(50,100), parameters = c(0.7, 0.3), nbTrials = 5, type = "binom")
#' dataGenerator_1D(chpts = c(50,100), parameters = c(0.4,0.7), nbSuccess = 10, type = "negbin")
#' dataGenerator_1D(chpts = c(50,70,120,200), parameter = c(0,3,-1,1), type = "gauss")
#' dataGenerator_1D(chpts = c(50,100,180), parameter = c(3,1,6), type = "variance")
dataGenerator_1D <- function(chpts = 100,
                             parameters = 0.5, # as 0 is not available for costs exp, poisson, geom, negbin
                             sdNoise = 1,
                             gamma = 1,
                             nbTrials = 10,
                             nbSuccess = 10,
                             type = "gauss")
{
  ############
  ### STOP ###
  ############

  if(!is.numeric(chpts)){stop('chpts values are not all numeric')}
  if(!all(chpts > 0)){stop('chpts values are not all positives')}
  if(is.unsorted(chpts, strictly = TRUE)){stop('chpts should be a strictly increasing vector of change-point positions (indices)')}

  if(!is.numeric(parameters)){stop('parameters values are not all numeric')}
  if(length(chpts) != length(parameters)){stop('chpts and parameters vectors are of different size')}

  allowed.types <- c("gauss", "exp", "poisson", "geom", "bern", "binom", "negbin", "variance")
  if(!type %in% allowed.types){stop('type must be one of: ', paste(allowed.types, collapse=", "))}

  ###################################
  ### Distribution specific stops ###
  ###################################

  if(type == "gauss")
  {
    if(length(sdNoise) > 1){stop('sdNoise should be length-1 vector')}
    if(!is.numeric(sdNoise)){stop('sdNoise value is not numeric')}
    if(sdNoise < 0){stop('sdNoise cannot be negative')}
    if(!is.numeric(gamma)){stop('gamma values are not all numeric')}
    if(any(gamma > 1 | gamma <= 0)){stop('gamma is not between 0 and 1 (0 excluded)')}
    if((length(gamma) != length(chpts)) && (length(gamma) != 1)){stop('vector gamma can be of length 1 or of the length of chpts')}
  }

  if(type == "binom")
  {
    if(length(nbTrials) > 1){stop('nbTrials should be length-1 vector')}
    if(!is.numeric(nbTrials)){stop('nbTrials value is not numeric')}
    if((nbTrials%%1 != 0) || (nbTrials <= 0)){stop('nbTrials cannot be non-positive or non-integer')}
  }

  if(type == "negbin")
  {
    if(length(nbSuccess) > 1){stop('nbSuccess should be length-1 vector')}
    if(!is.numeric(nbSuccess)){stop('nbSuccess value is not numeric')}
    if((nbSuccess%%1 != 0) || (nbSuccess <= 0)){stop('nbSuccess cannot be non-positive or non-integer')}
  }

  #############################
  ### parameter constraints ###
  #############################

  if(type == "variance"){if(min(parameters) < 0){stop('no negative parameters (variances) possible for Variance model')}}
  if(type == "exp"){if(min(parameters) <= 0){stop('no negative mean allowed for Poisson model')}}
  if(type == "poisson"){if(min(parameters) <= 0){stop('no negative mean allowed for Poisson model')}}

  if(type == "bern" || type == "binom")
    {if(any(parameters > 1) || any(parameters < 0)){stop('parameters should be probabilities between 0 and 1 (included)')}}

  if(type == "geom" || type == "negbin")
  {if(any(parameters > 1) || any(parameters <= 0)){stop('parameters should be probabilities between 0 and 1 (0 excluded)')}}

  ############  data generation   ############

  n <- chpts[length(chpts)]
  repetition <- c(chpts[1], diff(chpts))
  mu <- rep(parameters, repetition)

  if(type == "gauss" && all(gamma == 1)){y <- rnorm(n, mean = mu, sd = sdNoise)}

  if(type == "variance"){y <- rnorm(n, mean = 0, sd = mu)}

  if(type == "exp"){y <- rexp(n = n, rate = mu)}
  if(type == "poisson"){y <- rpois(n = n, lambda = mu)}
  if(type == "geom"){y <- rgeom(n = n, prob = mu) + 1} ### number of Bernoulli trials needed to get one success

  if(type == "bern"){y <- rbinom(n = n, size = 1, prob = mu)}
  if(type == "binom"){y <- rbinom(n = n, size = nbTrials, prob = mu)}
  if(type == "negbin"){y <- rnbinom(n = n, size = nbSuccess, prob = mu)}

  if(type == "gauss" && all(gamma < 1))
  {
    if(all(gamma < 1))
    {
      if(length(gamma) == 1){gamma <- rep(gamma, length(chpts))} #if gamma = 1 value, repeat it for each segment
      decay <- NULL
      for(i in 1:length(repetition)){decay <- c(decay, cumprod(rep(gamma[i], repetition[i]))/gamma[i])}
      y <- rnorm(n, mean = mu*decay, sd = sdNoise)
    }
  }
  return(y)
}

