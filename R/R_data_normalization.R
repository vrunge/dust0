

#########################################
#############  sdDiff_R  ################
#########################################

#' sdDiff_R
#'
#' @description sdDiff_R is a function based on the difference operator (or difference order for HALL method) estimating the time-series standard deviation in Gaussian model
#' The estimation is efficient for time-series with constant standard deviation and multiple changes in mean.
#' Three estimators are available:
#' \itemize{
#'   \item \code{HALL} : the so-called HALL-estimator of order 3. For more details see: \emph{(1990) Asymptotically optimal difference-based estimation of variance in nonparametric regression. Authors: Hall, Peter and Kay, JW and Titterinton, DM. Biometrika, pages 521--528}
#'   \item \code{MAD} : the Median Absolute Deviation estimator computed on \code{diff(y)/sqrt(2)} with y the vector of data-points
#'   \item \code{SD} : the standard deviation estimator (function sd) computed on \code{diff(y)/sqrt(2)} with y the vector of data-points
#' }
#' @param y Univariate time series
#' @param method Three available methods: \code{"HALL"}, \code{"MAD"} and \code{"SD"}
#' @return a value equal to the estimated standard deviation
#' @examples
#' data <- dataGenerator_1D(chpts = c(100, 200, 300), parameters = c(0,1,0), sdNoise = 2)
#' sdDiff_R(data, method = "HALL")
#' sdDiff_R(data, method = "MAD")
#' sdDiff_R(data, method = "SD")
sdDiff_R <- function(y, method = "HALL")
{
  if(is.numeric(y) == FALSE || length(y) < 5){stop('y is not a numeric vector or length < 5 (the HALL method cannot be used)')}
  if(method == "HALL")
  {
    n <- length(y)
    wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
    mat <- wei %*% t(y)
    mat[2, -n] <- mat[2, -1]
    mat[3, -c(n-1, n)] <- mat[3, -c(1, 2)]
    mat[4, -c(n-2, n-1, n)] <- mat[4, -c(1, 2, 3)]
    return(sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3)))
  }
  if(method == "MAD")
  {
    return(mad(diff(y)/sqrt(2)))
  }
  if(method == "SD")
  {
    return(sd(diff(y)/sqrt(2)))
  }
  return(NULL)
}
