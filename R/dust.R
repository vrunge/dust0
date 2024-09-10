
#' Run the 1D Dust Change Point Detection
#'
#' This function performs change point detection on one-dimensional data using the DUST algorithm.
#'
#' @param data A numeric vector. The time series data on which change point detection is performed.
#' @param penalty A numeric value. The penalty applied for adding a new change point. By default, it is set to \code{2 * log(length(data))}.
#' @param model A character string. Specifies the model used for change point detection. Default is "gauss". Possible values could include "gauss", "poisson", "exp", "geom", "binom", "negbin".
#' @param method A character string. Specifies the method to be used in the detection. Default is "fastest". Different segmentation methods may be specified.
#' @param alpha A numeric value. For randomness level. Default is \code{1e-9}.
#' @param nbLoops An integer. The number of loops to run in the max dual optimization algorithm. Default is 10.
#'
#' @return A list containing the change points detected by the DUST algorithm.
#'
#' @examples
#' data <- rnorm(1000)
#' result <- dust.1D(data = data)
#' @export
dust.1D <- function(
    data = data
    , penalty = 2*log(length(data))
    , model = "gauss"
    , method = "fastest"
    , alpha = 1e-9
    , nbLoops = 10
)
{
  partitioner <- new(DUST_1D, model, method, alpha, nbLoops)
  return(partitioner$quick_raw(data, penalty))
}
