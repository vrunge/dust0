

shift <- function(k){return(k+1)}

###
### mean from v_j to v_k included
###

eval_mean <- function(v, j, k){return((v[k+1]-v[j])/(k-j+1))}


###
###Variance for data-point y_j to y_k
###

eval_var <- function(cumData, cumData2, j, k)
{
  if(j == k){return(0)}
  return(eval_mean(cumData2, j, k) - eval_mean(cumData, j, k)^2)
}
