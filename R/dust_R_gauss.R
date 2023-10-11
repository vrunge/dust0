


dust_R_gauss <- function(data, penalty)
{
  #########
  ###
  ### DATA preprocessing
  ###
  n <- length(data)
  cumData <- cumsum(c(0, data))
  cumData2 <- cumsum(c(0, data^2))

  #########
  ###
  ### INITIALIZATION
  ###
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -penalty

  nb <- rep(0, n) #nb element for minimization in DP
  indexSet <- NULL


  #########
  ###
  ### update rule Dynamic Programming
  ###
  for(t in 1:n) # at t, transform Q_{t-1} into Q_{t}
  {
    indexSet <- c(indexSet, t) #add new test point

    min_temp <- Inf

    for(k in indexSet)
    {
      eval <- costQ[shift(k-1)] + (t-k+1)*eval_var(cumData, cumData2, k, t) + penalty
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t)] <- min_temp
    cp[shift(t)] <- index - 1

    #########
    ###
    ### PRUNING STEP PELT
    ###
    nonpruned <- NULL
    for(k in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inequality-based)
      if(costQ[shift(t)] > costQ[shift(k-1)] + (t-k+1)*eval_var(cumData, cumData2, k, t))
      {
        nonpruned <- c(nonpruned, k)
      }
    }
    indexSet <- nonpruned

    #########
    ###
    ### PRUNING STEP DUST
    ###
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    nonpruned <- NULL
    for(k in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inequality-based)
      if(costQ[shift(t)] > costQ[shift(k-1)] + (t-k+1)*eval_var(cumData, cumData2, k, t))
      {
        nonpruned <- c(nonpruned, k)
      }
    }
    indexSet <- nonpruned
    # TO DO
    # TO DO
    # TO DO
    nb[t] <- length(indexSet) ### count number of elements
  }

  #########
  ###
  ### backtracking step
  ###
  changepoints <- n # vector of change-point to build
  current <- n

  while(changepoints[1] > 0)
  {
    pointval <- cp[shift(current)] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  return(list(changepoints = changepoints[-1],  nb = nb, costQ = costQ))
}



