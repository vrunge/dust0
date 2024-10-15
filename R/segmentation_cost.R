
#' Compute the Total Segmentation Cost in One Dimension with Various Models
#'
#' This function calculates the total segmentation cost for one-dimensional data using various models such as Gaussian, Poisson, Exponential, and others.
#' It computes the segmentation cost for each segment defined by change-points (`chpts`) and sums these costs to provide a total segmentation cost.
#'
#' @param data A numeric vector representing the one-dimensional data to be segmented.
#' @param chpts A numeric vector of change-points. These indices represent the locations where the data is segmented. The vector should contain indices before any adjustments.
#' @param model A character string specifying the cost model to be used. Supported models include `"gauss"`, `"poisson"`, `"exp"`, `"geom"`, `"bern"`, `"binom"`, `"negbin"`, and `"variance"`. The default model is `"gauss"`.
#'
#' @return A numeric value representing the total segmentation cost.
#'
#' @details
#' For each segment defined by two consecutive change-points, the function applies \code{\link{Cost_1D}}, which calculates the cost of the segment based on the provided model. Supported models include:
#'
#' - `"gauss"`: Gaussian model, computes the negative log-likelihood under the Gaussian distribution.
#' - `"poisson"`: Poisson model.
#' - `"exp"`: Exponential model.
#' - `"geom"`: Geometric model.
#' - `"bern"`: Bernoulli model.
#' - `"binom"`: Binomial model.
#' - `"negbin"`: Negative Binomial model.
#' - `"variance"`: Variance-based cost model.
#'
#' If the model is `"gauss"`, the function adds an additional term `sum(data^2)/2` to the total cost to ensure that the cost is zero when the data is perfectly segmented with no noise.
#'
#' @examples
#' data <- rnorm(100)  # Generate some random Gaussian data
#' chpts <- c(25, 50, 75, 100)  # Example change-points
#' segmentation_Cost_1D(data, chpts, model = "gauss")
#'
#' @seealso \code{\link{Cost_1D}}
#'
#' @export
segmentation_Cost_1D <- function(data, chpts, model = "gauss")
{
  ### we add 0 and thus move all the indices
  chpts <- c(0, chpts) + 1
  K <- length(chpts) ### K >= 2
  S <- c(0, cumsum(data)) #WARNING: statistic is here alwoays the identity

  totalCost <- 0
  for (i in 2:K)
  {
    totalCost <- totalCost + Cost_1D(S, chpts[i-1], chpts[i], model)  # Apply the cost function
  }

  ### we add the sum of square in case of the Gaussian cast
  ### to get a 0 cost value in case of a no-noise perfectly well segmented data
  if(model == "gauss")
  {
    totalCost <- totalCost + sum(data^2)/2
  }
  return(totalCost)
}



#' Compute the Cost for a Single Segment Based on a Specified Model
#'
#' This function computes the cost of a single segment of data, defined by indices `a` and `b`,
#' using a model specified in the `model` parameter.
#'
#' @param S A cumulative sum statistic for the data.
#' @param a An integer representing the end index of the previous segment
#' @param b An integer representing the end index of the segment.
#' @param model A character string specifying the model to be used. Supported models include
#' `"gauss"`, `"poisson"`, `"exp"`, `"geom"`, `"bern"`, `"binom"`, `"negbin"`, and `"variance"`.
#'
#' @return A numeric value representing the cost of the segment.
#'
#' @details
#' The function supports several models, including:
#' - `"gauss"`: Gaussian model, which computes the negative log-likelihood for a Gaussian distribution.
#' - `"poisson"`: Poisson model
#' - `"exp"`: Exponential model
#' - `"geom"`: Geometric model
#' - `"bern"`: Bernoulli model
#' - `"binom"`: Binomial model
#' - `"negbin"`: Negative Binomial model
#' - `"variance"`: Variance-based cost model
#'
#' @export
Cost_1D <- function(S, a, b, model)
{

  cost_value <- 0
  diff <- S[b] - S[a]
  delta <- b - a

  if(model == "gauss")
  {
    cost_value <- - 0.5 * diff^2 / delta;
  }
  if (model == "poisson")
  {
    if (diff != 0.0)
    {
      cost_value <- diff * (1.0 - log(diff / delta))
    }
  }
  if (model == "exp")
  {
    cost_value <- delta * (1 + log(diff / delta))
  }
  if (model == "geom")
  {
    ratio <- diff / delta
    if(ratio != 1)
    {
      cost_value <- delta * log(ratio - 1) - diff * log((ratio - 1) / ratio)
    }
  }
  if (model == "bern")
  {
    ratio <- diff / delta
    if(ratio != 0 && ratio != 1)
    {
      cost_value <- - delta * (ratio * log(ratio) + (1 - ratio) * log(1 - ratio))
    }
  }
  if (model == "binom")
  {
    ratio <- diff / delta
    if(ratio != 0 && ratio != 1)
    {
      cost_value <- - delta * (ratio * log(ratio) + (1 - ratio) * log(1 - ratio))
    }
  }
  if (model == "negbin")
  {
    ratio <- diff / delta
    if(ratio != 0)
    {
      cost_value <- delta * log(1 + ratio) - diff * log(ratio / (1 + ratio))
    }
  }
  if (model == "variance")
  {
    cost_value <- 0.5 * delta * (1.0 + log(diff / delta))
  }
  return(cost_value)
}


