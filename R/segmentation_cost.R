


segmentation_Cost <- function(data, chpts, model = "gauss")
{
  chpts <- c(0, chpts)

  totalCost <- 0
  for (i in 2:length(chpts))
  {
    totalCost <- totalCost + Cost(model = model, chpts[i], chpts[i-1])  # Apply the cost function
  }

  return(totalCost)
}
