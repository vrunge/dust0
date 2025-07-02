
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


###########################################
#############  evalPrimal  ################
###########################################

#' evalPrimal
#'
#' @description Evaluation of the primal function
#' @param nat_theta the value of the natural parameter theta
#' @param A the nonlinear function
#' @param data the data to use in sum
#' @param const the constant term
evalPrimal <- function(nat_theta, A, data, const)
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

#' DD
#'
#' @description denominator value
#' @param mu the mu value
DD <- function(mu){return(1 - mu)}

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

  res1 <- (t-s1)* DD(mu) * (A(B(Ratio)) - Ratio* Bratio2)
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


##################################################
#############  evalDual_meanVar2 ################
##################################################

#' evalDual_meanVar2
#'
#' @description Evaluate the D(x) function as defined.
#' @param x The evaluation point
#' @param S Cumulative sum vector of data
#' @param S2 Cumulative sum vector of squared data
#' @param r First index
#' @param s Second index
#' @param t Final index
#' @param Qr The value of the global cost at r
#' @param Qs The value of the global cost at s
#' @param Qt The value of the global cost at t
evalDual_meanVar2 <- function(x, S, S2, r, s, t, Qr, Qs, Qt)
{
  # Segment lengths
  len_st <- t - s
  len_rs <- s - r

  # Mean of squared values
  mean_S2_st <- (S2[t] - S2[s]) / len_st
  delta_mean_S2_rst <- mean_S2_st - ((S2[s] - S2[r]) / len_rs)

  # Mean values
  mean_S_st <- (S[t] - S[s]) / len_st
  delta_mean_S_rst <-  mean_S_st - ((S[s] - S[r]) / len_rs)

  # Composite terms
  m2 <- mean_S2_st + x * delta_mean_S2_rst
  m1 <- mean_S_st + x * delta_mean_S_rst

  variance_term <- m2 - m1^2

  if (variance_term <= 0) {return(-Inf)}
  if (Qr == Inf) {return(-Inf)}
  mean_Q_st  <- (Qt - Qs) / len_st
  mean_Q_rst  <- mean_Q_st - ((Qs - Qr) / len_rs)

  res1 <- 0.5 * (1 + log(variance_term))
  res2 <- mean_Q_st + x * mean_Q_rst

  return(res1 - res2)
}


dualVALUE_meanVar2D <- function(S, S2,
                                r1, r2, s, t,
                                Qr1, Qr2, Qs, Qt)
{
  # Segment lengths
  len_st <- t - s
  len_r1s <- s - r1
  len_r2s <- s - r2

  # Mean of squared values
  A <- (S2[t] - S2[s]) / len_st
  D <- A - ((S2[s] - S2[r1]) / len_r1s)
  Dbar <- A - ((S2[s] - S2[r2]) / len_r2s)

  # Mean values
  a <- (S[t] - S[s]) / len_st
  d <-  a - ((S[s] - S[r1]) / len_r1s)
  dbar <-  a - ((S[s] - S[r2]) / len_r2s)

  delta  <- (Qt - Qs) / len_st
  DELTA  <- delta - ((Qs - Qr1) / len_r1s)
  DELTAbar  <- delta - ((Qs - Qr2) / len_r2s)

  #y = alpha + x beta
  ratio <- (Dbar*DELTA - D*DELTAbar)/(2*dbar*DELTA-2*d*DELTAbar)

  alpha <- A  - ratio^2 +  (Dbar/dbar) * (ratio - a)
  beta <- D - Dbar*d/dbar
  gamma0 <- delta  + (DELTAbar/dbar) * (ratio - a)
  gamma1 <- DELTA  - DELTAbar*d/dbar
  #print("res")
  #print(c(alpha, beta, gamma0, gamma1))
  #print(beta/(2*gamma1))
  if(Qr1 == Inf){gamma0 <- NA}
  if(Qr2 == Inf){gamma0 <- NA}
  if(Qs == Inf){gamma0 <- NA}
  u <- -d/dbar
  v <- (1/dbar)*(ratio - a)
  return(c(alpha, beta, gamma0, gamma1, u, v)) # y = ux * v
}



#' compute_xstar
#'
#' @description Computes the critical point \eqn{x^\star}
#' This function uses empirical variances and mean differences between segments \eqn{[r,s]} and \eqn{[s,t]}
#'
#' @param S Numeric vector. Cumulative sum of data (e.g., \code{cumsum(y)} with a prepended 0).
#' @param S2 Numeric vector. Cumulative sum of squares (e.g., \code{cumsum(y^2)} with a prepended 0).
#' @param r First index
#' @param s Second index
#' @param t Final index
#' @param Qr The value of the global cost at r
#' @param Qs The value of the global cost at s
#' @param Qt The value of the global cost at t
#'
#' @return The optimal value \eqn{x^\star} as a numeric scalar. Returns 0 if \code{Qr == Inf}.
#' @export
compute_xstar <- function(S, S2, r, s, t, Qr, Qs, Qt)
{
  if(Qr == Inf){return(0)}
  len_st <- t - s
  len_sr <- s - r

  # Empirical variances
  mean_st <- (S[t] - S[s]) / len_st
  var_st <- (S2[t] - S2[s]) / len_st - mean_st^2

  mean_sr <- (S[s] - S[r]) / len_sr
  var_sr <- (S2[s] - S2[r]) / len_sr - mean_sr^2

  # Delta mean slope
  delta_mean <- mean_sr - mean_st
  delta_mean_sq <- delta_mean^2

  # x0 and x1
  x0 <- 0.5 * ((var_st - var_sr) / delta_mean_sq - 1)
  x1 <- x0^2 + var_st / delta_mean_sq

  # x2 (delta Q)

    q_st <- (Qt - Qs) / len_st
    q_sr <- (Qs - Qr) / len_sr
    x2 <- q_st - q_sr

  # Compute x_star
  inv_2x2 <- 1 / (2 * x2)
  root_term <- sqrt(x1 + inv_2x2^2)
  x_star <- max(0, x0 + inv_2x2 - sign(x2) * root_term)
  return(x_star)
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



















