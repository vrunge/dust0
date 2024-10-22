
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
#' @param type the model: \code{"gauss"}, \code{"exp"}, \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}
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
#' dataGenerator_1D(chpts = c(50,70, 120, 200), parameter = c(0,3,-1,1), type = "gauss")
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



########################################################
#############    dataGenerator_MultiD   ################
########################################################

#' dataGenerator_MultiD
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
#' dataGenerator_MultiD(chpts = c(50,100),
#'                      parameters = data.frame(ts1 = c(10,10), ts2 = c(5,20)),
#'                      gamma = data.frame(ts1 = c(0.9,0.8), ts2 = c(0.8,0.9)),
#'                      sdNoise = 0.2, type = "gauss")
#' dataGenerator_MultiD(chpts = c(50,100,150),
#'                      parameters = data.frame(ts1 = c(4,0,2), ts2 = c(5,2,-2),ts3 = c(0,1,0)),
#'                      sdNoise = 0.4, type = "gauss")
#' dataGenerator_MultiD(chpts = c(50,100,150),
#'                      parameters = data.frame(ts1 = c(4,10,2), ts2 = c(5,2,2),ts3 = c(10,4,3)),
#'                      type = "poisson")
#' dataGenerator_MultiD(chpts = c(50,100,150),
#'      parameters = data.frame(ts1 = c(0.4,0.3,0.5), ts2 = c(0.5,0.6,0.2),ts3 = c(0.1,0.4,0.6)),
#'      nbTrials = c(30,10,100), type = "binom")
dataGenerator_MultiD <- function(chpts = 100,
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





########
########
########
########
########
########
########
########


######################################################
###########     dataGenerator_meanVar   ##############
######################################################

#' dataGenerator_meanVar
#'
#' @description Generating data for detecting changes in mean and variance with Gaussian model
#' @param chpts a vector of increasing change-point indices
#' @param means vector of successive segment means
#' @param sds vector of successive segment standard deviation
#' @return a vector of (univariate) simulated data with changes in mean and variance
#' @examples
#' dataGenerator_meanVar(chpts = c(30,100,120), means = c(0,1,0), sds = c(1,1,2))
dataGenerator_meanVar <- function(chpts = 100,
                                  means = 0,
                                  sds = 1)
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(chpts)){stop('chpts values are not all numeric')}
  if(!all(chpts > 0)){stop('chpts values are not all positives')}
  if(is.unsorted(chpts, strictly = TRUE)){stop('chpts should be a strictly increasing vector of change-point positions (indices)')}

  if(!is.numeric(means)){stop('means values are not all numeric')}
  if(length(chpts) != length(means)){stop('chpts and means vectors are of different size')}
  if(length(chpts) != length(sds)){stop('chpts and sds vectors are of different size')}
  if(min(sds) < 0){stop('sds values not all non-negative')}


  ############  data generation  ############

  n <- chpts[length(chpts)]
  repetition <- c(chpts[1], diff(chpts))

  mu <- rep(means, repetition)
  sds <- rep(sds, repetition)
  y <- rnorm(n, mean = mu, sd = sds)
  return(y)
}



########################################################
#############     dataGenerator_Reg    #################
########################################################

#' dataGenerator_Reg
#'
#' @description Generating data for changes in simple regression model and simple logistic regression model
#' @param chpts a vector of increasing change-point indices
#' @param A vector of regression coefficients A in Ax+B simple regression model
#' @param B vector of regression coefficients B in Ax+B simple regression model
#' @param meansX vector of mean values for x values generated by a Gaussian model (if one value, we copy the value to match chpts length)
#' @param sdX vector of standard deviation values for x values generated by a Gaussian model  (if one value, we copy the value to match chpts length)
#' @param sdNoise standard deviation of the Gaussian noise
#' @param type \code{"simple"} for the simple regression, \code{"logistic"} for the logistic regression
#' @return a dataframe with time series x and y for change-points in regression of type \code{y = Ax + B + noise} or\code{P(Y = 1) = 1/(1+exp(-(Ax+B+noise)))}.  \code{x} and \code{noise} are generated by a random Gaussian model.
#' @examples
#' dataGenerator_Reg(chpts = c(40,90), A = c(2,-1),  B = c(-1,2), meansX = c(1,2), type = "simple")
#' dataGenerator_Reg(chpts = c(50,100), A = c(-1,1),  B = c(1,-1), type = "logistic")
dataGenerator_Reg <- function(chpts = 100,
                              A = 2,
                              B = -1,
                              meansX = 0,
                              sdX = 1,
                              sdNoise = 1,
                              type = "simple")
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(chpts)){stop('chpts values are not all numeric')}
  if(!all(chpts > 0)){stop('chpts values are not all positives')}
  if(is.unsorted(chpts, strictly = TRUE)){stop('chpts should be a strictly increasing vector of change-point positions (indices)')}

  if(length(chpts) != length(A)){stop('chpts and A vectors are of different size')}
  if(length(chpts) != length(B)){stop('chpts and B vectors are of different size')}
  if(min(sdX) < 0){stop('sdX cannot be negative')}

  if(length(meansX) == 1){meansX <- rep(meansX, length(chpts))}
  if(length(sdX) == 1){sdX <- rep(sdX, length(chpts))}
  if(length(chpts) != length(meansX)){stop('chpts and meansX vectors are of different size')}
  if(length(chpts) != length(sdX)){stop('chpts and sdX vectors are of different size')}

  ############  data generation  ############

  n <- chpts[length(chpts)]
  repetition <- c(chpts[1], diff(chpts))
  A_rep <- rep(A, repetition)
  B_rep <- rep(B, repetition)
  meansX_rep <- rep(meansX, repetition)
  sdX_rep <- rep(sdX, repetition)

  x <- rnorm(n, mean = meansX_rep, sd = sdX_rep)

  if(type == "simple")
  {
    y <- A_rep*x + B_rep + rnorm(n, mean = 0, sd = sdNoise)
  }
  if(type == "logistic")
  {
    z <- A_rep*x + B_rep + rnorm(n, mean = 0, sd = sdNoise)
    proba <- 1/(1 + exp(-z))
    y <- rbinom(n, 1, proba)
  }

  return(data.frame(x = x, y = y))
}


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



