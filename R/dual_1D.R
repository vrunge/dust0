



dual_1D <- function(mu, x, s1, s2, t, type = "gauss", OP = FALSE, penalty = 2*length(x))
{
  ### preprocessing
  n <- length(x)
  x <- statistic(x)
  S <- cumsum(x)

  # GOAL : comparing q_t^s1 against q_t^s2 for pruning s1 with dust method

  # loading the type specific functions
  A <- A(type = type)
  B <- B(type = type)

  # mu rescaling
  mu <- mu * mu_max(S, s1, s2, t, type)  #mu was between 0 and 1

  if(OP == FALSE)
  {
    res <- evalDual(mu, A, B, S, s1, s2, 0, 0) #TO DO should work for mu = vector!!!
  }

  if(OP == TRUE)
  {
    OPres <- OP_R(x, penalty, type = type)
    const <- OPres$costQ
    res <- evalDual(mu, A, B, S, s1, s2, const[s1], const[s2]) #TO DO should work for mu = vector!!!
  }

  return(res)
}
