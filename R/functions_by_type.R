

statistic <- function(type)
{
  if(type == "gauss"){statistic <- function(x) x}
  if(type == "bern"){statistic <- function(x) x}
  if(type == "poisson"){statistic <- function(x) x}
  if(type == "exp"){statistic <- function(x) x}
  return(statistic)
}


###
### A
###

A <- function(type)
{
  if(type == "gauss"){A <- function(theta) theta^2/2}
  if(type == "bern"){A <- function(theta) log(1 + exp(theta))}
  if(type == "poisson"){A <- function(theta) exp(theta)}
  if(type == "exp"){A <- function(theta) -log(-theta)}
  return(A)
}

###
### B = (grad A)^(-1)
###

B <- function(type)
{
  if(type == "gauss"){B <- function(theta) theta}
  if(type == "bern"){B <- function(theta) log(1 + exp(theta))}
  if(type == "poisson"){B <- function(theta) exp(theta)}
  if(type == "exp"){B <- function(theta) -log(-theta)}
  return(B)
}

###
### mu_max
###

mu_max <- function(S, s1, s2, t, type) ### TO DO !!!
{
  if(type == "gauss"){res <- 1}
  if(type == "bern"){res <- 1}
  if(type == "poisson"){res <- 1}
  if(type == "exp"){res <- 1}
  return(res)
}


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

D <- function(mu, s1, s2, n){return((n - s1) - mu*(s1 - s2))}
R <- function(mu, S, s1, s2, n){return(((S[n] - S[s1]) - mu * (S[s1] - S[s2])) / ((n - s1) - mu * (s1 - s2)))}


evalDual <- function(mu, A, B, S, s1, s2, const1, const2)
{
  n <- length(S)
  res1 <- D(mu, s1, s2, n) * (A(B(R(mu, S, s1, s2, n))) - R(mu, S, s1, s2, n) * B(R(mu, S, s1, s2, n)))
  res2 <- const1 + mu*(const1 - const2)
  return(res1 + res2)
}


