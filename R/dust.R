
#' Run the 1D Dust Change Point Detection
#'
#' This function performs change point detection on one-dimensional data using the DUST algorithm.
#'
#' @param data A numeric vector. The time series data on which change point detection is performed.
#' @param penalty A numeric value. The penalty applied for adding a new change point. By default, it is set to \code{2 * log(length(data))}.
#' @param model A character string. Specifies the model used for change point detection. Default is "gauss". Possible values could include "gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance".
#' @param method A character string specifying the method used to handle indices and pruning tests in the algorithm. The default is \code{fastest}, which automatically selects the quickest method for the chosen model. Other available methods are:
#' \itemize{
#'   \item \code{"randIndex_Eval0"} to \code{"randIndex_Eval5"}: Random index-based methods with different dual maximization algorithm (0 through 5).
#'   \item \code{"detIndex_Eval0"} to \code{"detIndex_Eval5"}: Deterministic index-based methods  with different dual maximization algorithm (0 through 5).
#' }
#' Here are the current available algorithms (\code{Eval4} is often the most efficient one)
#' \itemize{
#'   \item \code{"Eval0"}: random evaluation of the dual (with uniform distribution)
#'   \item \code{"Eval1"}: max value with closed formula (gauss model only), otherwise no pruning performed and we get the (slow) OP algorithm
#'   \item \code{"Eval2"}: golden-section search.
#'   \item \code{"Eval3"}: binary search. At each step, we evaluate the tangent line to the current point at its max to stop the search at early step (when possible)
#'   \item \code{"Eval4"}: auasi-Newton method with armijo condition
#'   \item \code{"Eval5"}: PELT rule
#'   \item \code{"Eval6"}: OP rule
#' }
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





#' Run the MD Dust Change Point Detection
#'
#' This function performs change point detection on one-dimensional data using the DUST algorithm.
#'
#' @param data A numeric matrix. The time series data on which change point detection is performed. Each row should be a different dimension.
#' @param penalty A numeric value. The penalty applied for adding a new change point. By default, it is set to \code{2 nrow(data) log(length(data))}.
#' @param constraints An integer. The maximum number of constraints to be considered in the DUST pruning test.
#' @param model A character string. Specifies the model used for change point detection. Default is "gauss". Possible values could include "gauss", "poisson", "exp", "geom", "binom", "negbin", "variance".
#' @param method A character string. Specifies the method to be used in the detection. Default is "fastest". Different segmentation methods may be specified.
#' @param alpha A numeric value. For randomness level. Default is \code{1e-9}.
#' @param nbLoops An integer. The number of loops to run in the max dual optimization algorithm. Default is 10.
#'
#' @return A list containing the change points detected by the DUST algorithm.
#'
#' @examples
#' data <- matrix(rnorm(1000 * 2), nrow = 2, ncol = 1000)
#' result <- dust.MD(data = data)
#' @export
dust.MD <- function(
    data = data
    , penalty = 2*nrow(data)*log(length(data))
    , constraints = nrow(data)
    , model = "gauss"
    , method = "fastest"
    , alpha = 1e-9
    , nbLoops = 10
)
{
  partitioner <- new(DUST_MD, model, method, alpha, nbLoops)
  return(partitioner$quick_raw(data, penalty, constraints))
}




