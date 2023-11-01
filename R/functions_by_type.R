

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

##########################################################################################
###
### A
###

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

##########################################################################################
###
### B = (grad A)^(-1)
###

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

##########################################################################################
###
### mu_max (work for 1D case only)
###

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


##########################################################################################
###
### EVAL FUNCTION
###

eval <- function(nat_theta, A, data, const)
{
  ### data transformed by function statistic (often = identity)
  ### nat_theta = natural parameter (not theta, but its transformation)
  ###             it range : definition set of A
  return(length(data)*A(nat_theta) - sum(data)*nat_theta + const)
}

##########################################################################################
###
### MINIMAL COST
###

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


##########################################################################################
###
### denominator and ratio
###

D <- function(mu){return(1 - mu)}
R <- function(mu, S, s1, s2, t)
{
  Mt <- (S[t] - S[s1])/(t-s1)
  Ms <- (S[s1] - S[s2])/(s1-s2)
  if(Mt == Ms){return(Mt)}
  return((Mt - mu * Ms) / (1 - mu))
}


##########################################################################################
###
### evalDual
###
###
### DANGER : don't use it at mu_max value
###

evalDual <- function(mu, A, B, S, s1, s2, t, const1, const2)
{
  Ratio <- R(mu, S, s1, s2, t)

  Bratio2 <- B(Ratio)
  Bratio2[is.infinite(Bratio2)] <- 0 ### to solve Ratio * B(Ratio) = 0 * Inf

  res1 <- (t-s1)* D(mu) * (A(B(Ratio)) - Ratio* Bratio2)
  res2 <- const1 + mu * ((t-s1)/(s1-s2)) * (const1 - const2)

  return(res1 + res2)
}


