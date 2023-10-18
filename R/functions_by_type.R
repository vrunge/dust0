

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


###
### A
###

A <- function(type)
{
  if(type == "gauss"){A <- function(theta) theta^2/2}
  if(type == "exp"){A <- function(theta) -log(-theta)}
  if(type == "poisson"){A <- function(theta) exp(theta)}
  if(type == "geom"){A <- function(theta) -log(exp(-theta)-1)}
  if(type == "bern"){A <- function(theta) log(1 + exp(theta))}
  if(type == "binom"){A <- function(theta) log(1 + exp(theta))}
  if(type == "negbin"){A <- function(theta) -log(1 - exp(theta))}
  return(A)
}

###
### B = (grad A)^(-1)
###

B <- function(type)
{
  if(type == "gauss"){B <- function(theta) theta}
  if(type == "exp"){B <- function(theta) -1/theta}
  if(type == "poisson"){B <- function(theta) log(theta)}
  if(type == "geom"){B <- function(theta) theta <- log((theta-1)/theta)}
  if(type == "bern"){B <- function(theta){theta <- pmax(pmin(theta,1),0);log(theta/(1-theta))}}
  if(type == "binom"){B <- function(theta){theta <- pmax(pmin(theta,1),0); log(theta/(1-theta))}}
  if(type == "negbin"){B <- function(theta){theta <- pmax(theta,0); log(theta/(1+theta))}}
  return(B)
}

###
### mu_max (work for 1D and MultiD cases)
###

mu_max <- function(S, s1, s2, t, type)
{
  if(is.vector(S)){St <- S[t]; Ss1 <- S[s1]; Ss2 <- S[s2]}
  if(is.matrix(S)){St <- S[,t]; Ss1 <- S[,s1]; Ss2 <- S[,s2]}

  if(type == "gauss"){res <- 1}
  if(type == "exp"){res <- min(1, (s1-s2)/(t-s1)*(St - Ss1) / (Ss1 - Ss2), na.rm = TRUE)}
  if(type == "poisson"){res <- min(1, (s1-s2)/(t-s1)*(St - Ss1) / (Ss1 - Ss2), na.rm = TRUE)}
  if(type == "geom"){res <- min(1, (s1-s2)/(t-s1)*(St - Ss1 - (t-s1)) / (Ss1 - Ss2 - (s1-s2)), na.rm = TRUE)}
  if(type == "bern"){res <- (s1-s2)/(t-s1) * min((St - Ss1) / (Ss1 - Ss2),
                                                 (St - Ss1 - (t-s1)) / (Ss1 - Ss2 - (s1-s2)), na.rm = TRUE)}
  if(type == "binom"){res <- (s1-s2)/(t-s1) * min((St - Ss1) / (Ss1 - Ss2),
                                                 (St - Ss1 - (t-s1)) / (Ss1 - Ss2 - (s1-s2)), na.rm = TRUE)}
  if(type == "negbin"){res <- min(1, (s1-s2)/(t-s1)*(St - Ss1) / (Ss1 - Ss2), na.rm = TRUE)}

  res <- max(res,  0)* (t-s1)/(s1-s2) ### NORMALIZATION (because not done in D and R), max to avoid case res = -Inf (bern, binom)

  return(res)
}

##########################################################################################
##########################################################################################
##########################################################################################

###
### EVAL FUNCTION
###

eval <- function(x, A, data, const)
{
 return(length(data)*A(x) - sum(data)*x + const)
}


###
### EVAL DUAL FUNCTION
###

D <- function(mu, s1, s2, t){return((t - s1) - mu*(s1 - s2))}
R <- function(mu, S, s1, s2, t){return(((S[t] - S[s1]) - mu * (S[s1] - S[s2])) / ((t - s1) - mu * (s1 - s2)))}


evalDual <- function(mu, A, B, S, s1, s2, t, const1, const2)
{
  Ratio <- R(mu, S, s1, s2, t)
  Bratio <- B(Ratio)
  Bratio2 <- Bratio
  Bratio2[is.infinite(Bratio2)] <- 0
  res1 <- D(mu, s1, s2, t) * (A(Bratio) - Ratio* Bratio2)
  res2 <- const1 + mu*(const1 - const2)
  return(res1 + res2)
}

###
### MINIMAL COST
###

min_cost <- function(A, B, S, s, t, const)
{
  data <- S[t] - S[s]
  point <- data/(t-s)
  value1 <- B(point)
  if(data == 0){value2 <- 0}else{value2 <- value1} #to avoid data*value = 0 * Inf
  res <- (t-s)*A(value1) - data*value2 + const
  if(is.nan(res)){res <- const} # cases geom, bern, binom with point = 1
  return(res)
}



