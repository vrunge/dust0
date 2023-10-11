

####################################################
#############    dataGenerator_1D   ################
####################################################


#' dataGenerator_1D
#'
#' @description Generating data for multiple change-point detection in univariate (uni-parametric) time series
#' @param chpts a vector of increasing change-point indices
#' @param parameters vector of successive segment parameters
#' @param sdNoise standard deviation for the noise parameter (type gauss)
#' @param nbTrials number of trials (type binom)
#' @param nbFailures number of failures (type negbin)
#' @param type the model: gauss, poisson
#' @return a vector of simulated time series
#' @examples
#' myData <- dataGenerator_1D(chpts = c(50,100), parameter = c(0,1), type = "gauss")
dataGenerator_1D <- function(chpts = 100,
                          parameters = 0,
                          sdNoise = 1,
                          nbTrials = 10,
                          nbFailures = 10,
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


  ###################################
  ### Distribution specific stops ###
  ###################################

  # TO DO
  if(sdNoise < 0){stop('sdNoise cannot be negative')}
  if(type == "poisson"){if(min(parameters) <= 0){stop('no negative mean allowed for Poisson model')}}

  ############
  ############
  ############

  n <- chpts[length(chpts)]
  repetition <- c(chpts[1], diff(chpts))
  mu <- rep(parameters, repetition)

  if(type == "gauss"){y <- rnorm(n, mean = mu, sd = sdNoise)}

  if(type == "exp"){y <- rexp(n = n, rate = mu)}
  if(type == "poisson"){y <- rpois(n = n, lambda = mu)}

  if(type == "negbin"){y <- rnbinom(n = n, size = 1, prob = mu)} # TO DO
  if(type == "geom"){y <- rgeom(n = n, prob = mu)}

  if(type == "bern"){y <- rbinom(n = n, size = 1, prob = mu)}
  if(type == "binom"){y <- rbinom(n = n, size = nbTrials, prob = mu)}

  if(type == "chi2"){y <- rchisq(n = n, df = mu)}

  if(type == "laplace"){y <- NULL} #write my function
  if(type == "pareto"){y <- NULL} #write my function

  return(y)
}





########################################################
#############    dataGenerator_MultiD   ################
########################################################


#' dataGenerator_MultiD
#'
#' @description Generating data for multiple change-point detection in multivariate time series binding together 1D time-series from dataGenerator1D
#' @param chpts a vector of increasing change-point indices
#' @param parameters matrix of successive segment parameters
#' @param sdNoise standard deviation for the noise parameter (type gauss)
#' @param nbTrials number of trials (type binom)
#' @param nbFailures number of failures (type negbin)
#' @param type the model: gauss, poisson
#' @return a vector of simulated time series
#' @examples
#' myData <- dataGenerator_MultiD(chpts = c(50,100), parameter = c(0,1), type = "gauss")
dataGenerator_MultiD <- function(chpts = 100,
                            parameters = 0,
                            sdNoise = 1,
                            nbTrials = 10,
                            nbFailures = 10,
                            type = "gauss")
{
  # call p times the function dataGenerator1D
}



######################################################
#############     dataGenerator_MV   #################
######################################################

#' dataGenerator_MV
#'
#' @description Generating data for detecting changes in mean and variance with Gaussian model
#' @param chpts a vector of increasing change-point indices
#' @param means vector of successive segment means
#' @param sds vector of successive segment standard deviation
#' @return a vector of (univariate) simulated data with changes in mean and variance
#' @examples
#' myData <- dataGenerator_MV(chpts = c(30,100,120), means = c(0,1,0), sds = c(1,1,2))
dataGenerator_MV <- function(chpts = 100,
                             means = 0,
                             sds = 1)
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(chpts)){stop('chpts values are not all numeric')}
  if(is.unsorted(chpts)){stop('chpts should be an increasing vector of change-point positions (indices)')}
  if(length(unique(chpts)) < length(chpts)){stop('chpts is not a strictly increasing sequence')}
  if(!is.numeric(means)){stop('means values are not all numeric')}
  if(length(chpts) != length(means)){stop('chpts and means vectors are of different size')}
  if(length(chpts) != length(sds)){stop('chpts and sds vectors are of different size')}
  if(!all(sds >= 0)){stop('sds not all non-negative')}

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
#' @description Generating data for changes in bivariate independent time series
#' @param chpts a vector of increasing change-point indices
#' @param A vector of regression coefficients A in A*x+B simple regression model
#' @param B vector of regression coefficients B in A*x+B simple regression model
#' @param meansX vector of mean values for x values generated by a Gaussian model
#' @param sdX vector of standard deviation values for x values generated by a Gaussian model
#' @param sdNoise standard deviation of the Gaussian noise
#' @return a dataframe with time series x and y for change-points in regression of type y = A*x + B + noise
#' @examples
#' myData <- dataGenerator_Reg(chpts = c(40,90), A = c(2,-1),  B = c(-1,2), meansX = c(1,2))
dataGenerator_Reg <- function(chpts = 100,
                              A = 2,
                              B = -1,
                              meansX = 0,
                              sdX = 1,
                              sdNoise = 1)
{
  ############
  ### STOP ###
  ############
  if(!is.numeric(chpts)){stop('data values are not all numeric')}
  if(is.unsorted(chpts)){stop('chpts should be an increasing vector of change-point positions (indices)')}
  if(length(unique(chpts)) < length(chpts)){stop('chpts is not a strictly increasing sequence')}
  if(length(chpts) != length(A)){stop('chpts and A vectors are of different size')}
  if(length(chpts) != length(B)){stop('chpts and B vectors are of different size')}
  if(min(sdX) < 0){stop('sdX cannot be negative')}

  if(length(meansX) == 1){meansX <- rep(meansX, length(chpts))}
  if(length(sdX) == 1){sdX <- rep(sdX, length(chpts))}
  if(length(chpts) != length(meansX)){stop('chpts and meansX vectors are of different size')}
  if(length(chpts) != length(sdX)){stop('chpts and sdX vectors are of different size')}

  n <- chpts[length(chpts)]
  repetition <- c(chpts[1], diff(chpts))
  A_rep <- rep(A, repetition)
  B_rep <- rep(B, repetition)
  meansX_rep <- rep(meansX, repetition)
  sdX_rep <- rep(sdX, repetition)

  x <- rnorm(n, mean = meansX_rep, sd = sdX_rep)
  y <- A_rep*x + B_rep + rnorm(n, mean = 0, sd = sdNoise)
  return(data.frame(x = x, y = y))
}


