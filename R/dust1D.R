
#'  Multiple Change Point Detection for 1D data with DUST algorithm
#'
#' This function performs change-point detection on univariate time series using the DUST algorithm
#'
#' @param data A numeric vector. The time series on which change-point detection is performed. There is no data copy and no possibility to append new data to complete the analysis (to do so, see \code{\link{dust.object.1D}})
#' @param penalty A positive numeric value. The penalty applied for adding a new change point. By default, it is set to \code{2 sdDiff(data) log(length(data))}.
#' @param model A character string. Specifies the model used for change-point detection. Default is "gauss". Possible values could include "gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance".
#' @param method A character string specifying the method used to handle indices and pruning tests in the algorithm. The default is \code{detIndex_Eval4}, which automatically selects the quickest method for the chosen model. Other available methods are:
#' \itemize{
#'   \item \code{"randIndex_Eval0"} to \code{"randIndex_Eval6"}: Random index-based methods with different dual maximization algorithm (0 through 6).
#'   \item \code{"detIndex_Eval0"} to \code{"detIndex_Eval6"}: Deterministic index-based methods  with different dual maximization algorithm (0 through 6).
#' }
#' Here are the current available algorithms for evaluation the maximum of the dual function (\code{Eval4} is often the most efficient one). In all cases, the samllest remaining index is tested at each iteration with PELT.
#' \itemize{
#'   \item \code{"Eval0"}: random evaluation of the dual (with uniform distribution)
#'   \item \code{"Eval1"}: max value with closed formula (gauss model only), otherwise no pruning performed and we get the (slow) OP algorithm
#'   \item \code{"Eval2"}: golden-section search.
#'   \item \code{"Eval3"}: binary search. At each step, we evaluate the tangent line to the current point at its max to stop the search at early step (when possible)
#'   \item \code{"Eval4"}: quasi-Newton method with armijo condition
#'   \item \code{"Eval5"}: PELT rule
#'   \item \code{"Eval6"}: OP rule
#' }
#' @param nbLoops An integer. The number of loops to run in the max dual optimization algorithm. Default is 10.
#'
#' @return A list containing the change points detected by the DUST algorithm and other info.
#'
#' @note The input data should be first normalized by function \code{data_normalization_1D} to use the default penalty in Gaussian model, instead of value `2*sdDiff(data)^2*log(length(data))`
#'
#' @seealso  \code{\link{dataGenerator_1D}} for easy 1D data generation with change points and different data types/models
#'   \code{\link{data_normalization_1D}} for normalizing data before using \code{dust.1D}
#'   \code{\link{dust.object.1D}} A function similar to \code{dust.1D} but explicitly utilizes a class-based structure. This approach copies data into an object, allowing incremental analysis by appending new data with \code{append_c} and  \code{update_partition}.
#' @examples
#' data <- rnorm(1000)
#' data <- data_normalization_1D(data)
#' result <- dust.1D(data = data)
#' @export
dust.1D <- function(
    data = data
    , penalty = 2*log(length(data))
    , model = "gauss"
    , method = "detIndex_Eval4"
    , nbLoops = 10
)
{
  partitioner <- new(DUST_1D, model, method, nbLoops)
  partitioner$one_dust(data, penalty)
  return(c(partitioner$get_partition(), partitioner$get_info()))
}


