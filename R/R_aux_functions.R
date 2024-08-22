
######################################
#############  shift  ################
######################################

#' shift
#'
#' @description Shifting integer by +1 value
#' @param k an integer
shift <- function(k){return(k+1)}


#########################################################
#############  backtracking_changepoint  ################
#########################################################

#' backtracking_changepoint
#'
#' @description Shifting integer by +1 value
#' @param cp vector of changepoints of size n+1
#' @param n data length
backtracking_changepoint <- function(cp, n)
{
  changepoints <- n ##### vector of change-point to build
  current <- n

  while(changepoints[1] > 0)
  {
    pointval <- cp[current] ##### new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }

  changepoints <- changepoints[-1] ##### remove the first point equal to 0

  return(changepoints)
}


