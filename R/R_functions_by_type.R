
##########################################
#############  statistic  ################
##########################################

#' statistic
#'
#' @description statistic function (cf exponential familly)
#' @param type the type of cost
statistic <- function(type)
{
  if(type == "gauss"){statistic <- function(x) x}
  if(type == "exp"){statistic <- function(x) x}
  if(type == "poisson"){statistic <- function(x) x}
  if(type == "geom"){statistic <- function(x) x}
  if(type == "bern"){statistic <- function(x) x}
  if(type == "binom"){statistic <- function(x) x}
  if(type == "negbin"){statistic <- function(x) x}
  return(statistic)
}



##################################
#############  A  ################
##################################

#' A
#'
#' @description Non linear function A
#' @param type the type of cost
A <- function(type)
{
  if(type == "gauss"){A <- function(nat_theta) nat_theta^2/2}
  if(type == "exp"){A <- function(nat_theta){u <- pmin(nat_theta,0); -log(-u)}} #Inf in nat_theta = 0 (nat_theta <= 0)
  if(type == "poisson"){A <- function(nat_theta) exp(nat_theta)}
  if(type == "geom"){A <- function(nat_theta){u <- pmin(nat_theta,0); -log(exp(-u)-1)}} #Inf in nat_theta = 0 (nat_theta <= 0)
  if(type == "bern"){A <- function(nat_theta) log(1 + exp(nat_theta))}
  if(type == "binom"){A <- function(nat_theta) log(1 + exp(nat_theta))}
  if(type == "negbin"){A <- function(nat_theta){u <- pmin(nat_theta,0); -log(1-exp(u))}} #Inf in nat_theta = 0 (nat_theta <= 0)
  return(A)
}


##################################
#############  B  ################
##################################

#' B
#'
#' @description Non linear function B = (grad A)^(-1)
#' @param type the type of cost
B <- function(type)
{
  ###
  ### here range theta = range of the data
  ###
  if(type == "gauss"){B <- function(theta) theta}
  if(type == "exp"){B <- function(theta){u <- pmax(theta,0); -1/u}} #-Inf in theta = 0 (theta >= 0)
  if(type == "poisson"){B <- function(theta){u <- pmax(theta,0); log(u)}} #-Inf in theta = 0 (theta >= 0)
  if(type == "geom"){B <- function(theta){u <- pmax(theta,1); theta <- log((u-1)/u)}} #-Inf in theta = 1 (theta >= 1)
  if(type == "bern"){B <- function(theta){u <- pmax(pmin(theta,1),0); log(u/(1-u))}}
                                        #-Inf in theta = 0, +Inf in theta = 1 (0 <= theta <= 1)
  if(type == "binom"){B <- function(theta){u <- pmax(pmin(theta,1),0); log(u/(1-u))}}
                                        #-Inf in theta = 0, +Inf in theta = 1 (0 <= theta <= 1)
  if(type == "negbin"){B <- function(theta){u <- pmax(theta,0); log(u/(1+u))}} #-Inf in theta = 0 (theta >= 0)
  return(B)
}


#######################################
#############  mu_max  ################ (work for 1D case only)
#######################################

#' mu_max
#'
#' @description Getting the argmax dual value
#' @param S The cumsum S vector
#' @param s1 the first index
#' @param s2 the second index
#' @param t the time step
#' @param type the type of cost
mu_max <- function(S, s1, s2, t, type)
{
  if(s2 > s1){return(Inf)}

  Mt <- (S[t] - S[s1])/(t - s1)
  Ms <- (S[s1] - S[s2])/(s1 - s2)

  if(type == "gauss"){res <- 1}
  if(type == "exp"){res <-  min(1,  Mt / Ms, na.rm = TRUE)}
  if(type == "poisson"){res <- min(1, Mt / Ms, na.rm = TRUE)}
  if(type == "geom"){res <- min(1, (Mt - 1) / (Ms - 1), na.rm = TRUE)}
  if(type == "bern"){res <-  min(Mt / Ms, (1 - Mt) / (1 - Ms), na.rm = TRUE)}
  if(type == "binom"){res <-  min(Mt / Ms, (1 - Mt) / (1 - Ms), na.rm = TRUE)}
  if(type == "negbin"){res <- min(1, Mt / Ms, na.rm = TRUE)}

  res <- max(res,  0) ### to avoid negative result (if res = -10^-16 for some reason...)

  return(res)
}


#####################################
#############  eval  ################
#####################################

#' eval
#'
#' @description Evaluation of the primal function
#' @param nat_theta the value of the natural parameter theta
#' @param A the nonlinear function
#' @param data the data to use in sum
#' @param const the constant term
eval <- function(nat_theta, A, data, const)
{
  ### data transformed by function statistic (often = identity)
  ### nat_theta = natural parameter (not theta, but its transformation)
  ###             it range : definition set of A
  return(length(data)*A(nat_theta) - sum(data)*nat_theta + const)
}

#########################################
#############  min_cost  ################
#########################################

#' min_cost
#'
#' @description minimal cost value
#' @param A the nonlinear function
#' @param B the B derived from A function
#' @param S the cumsum data
#' @param s the start index
#' @param t the end index
#' @param const the constant term
min_cost <- function(A, B, S, s, t, const)
{
  data <- S[t] - S[s]
  my_mean <- data/(t-s)
  argmin <- B(my_mean)

  if(data == 0){argmin2 <- 0}else{argmin2 <- argmin} ### to solve data*value = 0 * Inf
  res <- (t-s)*A(argmin) - data*argmin2 + const

  ###
  ### DANGER DANGER DANGER
  ###
  if(is.nan(res)){res <- const} # cases geom, bern, binom with my_mean = 1 (verified many times...)
  ###
  ### END OF DANGER DANGER DANGE
  ###

  return(res)
}




#######################################################
#############   denominator and ratio  ################
#######################################################

#' D
#'
#' @description denominator value
#' @param mu the mu value
D <- function(mu){return(1 - mu)}

#' R
#'
#' @description ratio
#' @param mu the mu value
#' @param S The cumsum S vector
#' @param s1 the first index
#' @param s2 the second index
#' @param t the time step
R <- function(mu, S, s1, s2, t)
{
  Mt <- (S[t] - S[s1])/(t-s1)
  Ms <- (S[s1] - S[s2])/(s1-s2)
  if(Mt == Ms){return(Mt)}
  return((Mt - mu * Ms) / (1 - mu))
}




##########################################
#############  evalDual   ################ DANGER : don't use it at mu_max value
##########################################

#' evalDual
#'
#' @description evaluation of the dual function
#' @param mu the mu value, evaluation point
#' @param A the nonlinear function
#' @param B the B derived from A function
#' @param S The cumsum S vector
#' @param s1 the first index
#' @param s2 the second index
#' @param t the time step
#' @param const1 first constant term
#' @param const2 second constant term
evalDual <- function(mu, A, B, S, s1, s2, t, const1, const2)
{
  Ratio <- R(mu, S, s1, s2, t)

  Bratio2 <- B(Ratio)
  Bratio2[is.infinite(Bratio2)] <- 0 ### to solve Ratio * B(Ratio) = 0 * Inf

  res1 <- (t-s1)* D(mu) * (A(B(Ratio)) - Ratio* Bratio2)
  res2 <- const1 + mu * ((t-s1)/(s1-s2)) * (const1 - const2)

  return(res1 + res2)
}



##################################################
#############  min_cost_meanVar   ################
##################################################

#' min_cost_meanVar
#'
#' @description evaluation of the cost function at its minimum for mean and var change Gaussian problem
#' @param S The cumsum S vector
#' @param S2 The second cumsum vector for squared data
#' @param s the start index
#' @param t the end step
#' @param const constant term
min_cost_meanVar <- function(S, S2, s, t, const)
{
  if(s + 1 == t){return(Inf)}
  Ms <- (S[t] - S[s])/(t-s)
  Ms2 <- (S2[t] - S2[s])/(t-s)
  res <- ((t-s)/2)*(1 + log(Ms2 - Ms^2)) + const
  return(res)
}

##################################################
#############  evalDual_meanVar   ################
##################################################

#' evalDual_meanVar
#'
#' @description evaluation of the dual function
#' @param mu the mu value, evaluation point
#' @param S The cumsum S vector
#' @param S2 The second cumsum vector for squared data
#' @param s1 the first index
#' @param s2 the second index
#' @param t the time step
#' @param const1 first constant term
#' @param const2 second constant term
evalDual_meanVar <- function(mu, S, S2, s1, s2, t, const1, const2)
{
  l <- (t - s1) - mu * (s1 - s2)
  D2 <- S2[t] - S2[s1] - mu * (S2[s1] - S2[s2])
  D1 <- S[t] - S[s1] - mu *(S[s1] - S[s2])

  res1 <- (1/2) * l * (1 + log((D2/l) - (D1/l)^2))
  res2 <- const1 + mu * (const1 - const2)

  return(res1 + res2)
}


#####################################################
#############  min_cost_regression   ################
#####################################################

#' min_cost_regression
#'
#' @description evaluation of the min cost value
#' @param A A
#' @param B B
#' @param C C
#' @param D D
#' @param E E
#' @param f f
#' @param s the start index
#' @param t the end index
#' @param const  constant term
min_cost_regression <- function(A, B, C, D, E, f, s, t, const)
{
  if(s + 1 == t){return(Inf)}
  Adiff <- A[t] - A[s]
  Bdiff <- B[t] - B[s]
  Cdiff <- C[t] - C[s]
  Ddiff <- D[t] - D[s]
  Ediff <- E[t] - E[s]
  Fdiff <- f[t] - f[s]

  num <- 2 * Bdiff * Ddiff * Ediff - Adiff* Ediff^2 - Cdiff * Ddiff^2
  denom <- Adiff * Cdiff - Bdiff^2
  res <- num/denom + Fdiff + const
  return(res)
}

#####################################################
#############  evalDual_regression   ################
#####################################################

#' evalDual_regression
#'
#' @description evaluation of the dual function
#' @param mu mu
#' @param A A
#' @param B B
#' @param C C
#' @param D D
#' @param E E
#' @param f f
#' @param s1 the first index
#' @param s2 the second index
#' @param t the time step
#' @param const1 first constant term
#' @param const2 second constant term
evalDual_regression <- function(mu,A,B,C,D,E,f, s1, s2, t, const1, const2)
{
  Adiff <- A[t] - A[s1] - mu * (A[s1] - A[s2])
  Bdiff <- B[t] - B[s1] - mu * (B[s1] - B[s2])
  Cdiff <- C[t] - C[s1] - mu * (C[s1] - C[s2])
  Ddiff <- D[t] - D[s1] - mu * (D[s1] - D[s2])
  Ediff <- E[t] - E[s1] - mu * (E[s1] - E[s2])
  Fdiff <- f[t] - f[s1] - mu * (f[s1] - f[s2])

  num <- 2 * Bdiff * Ddiff * Ediff - Adiff* Ediff^2 - Cdiff * Ddiff^2
  denom <- Adiff * Cdiff - Bdiff^2
  res1 <- num/denom + Fdiff
  res2 <- const1 + mu * (const1 - const2)
  return(res1 + res2)
}


###############################################
#############  mu_max_2param   ################
###############################################

#' mu_max_2param
#'
#' @description evaluation of the mu max in 2 param case
#' @param S The cumsum S vector
#' @param S2 The second cumsum vector for squared data
#' @param s1 the first index
#' @param s2 the second index
#' @param t the time step
mu_max_2param <- function(S, S2, s1, s2, t)
{
  if(s2 > s1){return(Inf)}

  Mt <- (S[t] - S[s1])/(t - s1)
  Mt2 <- (S2[t] - S2[s1])/(t - s1)
  Ms <- (S[s1] - S[s2])/(s1 - s2)
  Ms2 <- (S2[s1] - S2[s2])/(s1 - s2)

  Va <- Mt2 - Mt^2
  Vb <- Ms2 - Ms^2
  eps <- (Mt - Ms)/sqrt(Va + Vb)
  R <- (t-s1)/(s1-s2)
  u <- (Va+Vb)*(1+eps^2)
  delta <- u^2 - 4 * Va * Vb
  res <- R*(u - sqrt(delta))/(2*Vb)
  return(res)
}



















