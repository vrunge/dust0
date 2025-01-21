
#' Run the meanVar Dust Change Point Detection
#'
#' This function performs change point detection for mean and variance in univariate time series using the DUST algorithm.
#' @param data A numeric matrix. The time series data on which change point detection is performed. Each row represents a time-series.
#' @param penalty A numeric value. The penalty applied for adding a new change point. By default, it is set to \code{4 x log(length(data))}.
#' @param nbLoops An integer. The number of loops to run in the max dual optimization algorithm. Default is 10.
#' @export
dust.meanVar <- function(
      data = data
    , penalty = 2*nrow(data)*log(length(data))
    , nbLoops = 10
)
{
  partitioner <- new(DUST_meanVar, "", nbLoops)
  return(partitioner$quick_raw(data, penalty))
}


