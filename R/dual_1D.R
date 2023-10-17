



dual_1D <- function(mu, x, s1, s2, t, type = "gauss", OP = FALSE, penalty = 2*length(x))
{
  ###
  ### preprocessing
  ###
  n <- length(x)
  stat <- statistic(type = type)
  S <- cumsum(stat(x))
  ###
  ### GOAL : comparing q_t^s1 against q_t^s2 for pruning s1 with dust method
  ### loading the type specific functions
  A <- A(type = type)
  B <- B(type = type)

  ###
  ### mu-rescaling (to be between the right bounds, with respect to the type)
  ###
  mu <- mu * mu_max(S, s1, s2, t, type)  #mu was between 0 and 1

  if(OP == FALSE)
  {
    res <- evalDual(mu, A, B, S, s1, s2, t, 0, 0)
  }

  if(OP == TRUE)
  {
    OPres <- OP_R(x, penalty = penalty, type = type)
    cost <- OPres$costQ
    res <- evalDual(mu, A, B, S, s1, s2, t, cost[s1] + penalty, cost[s2] + penalty)
  }

  return(res)
}
