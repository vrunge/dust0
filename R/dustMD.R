if (!exists("DUSTMODULEMD_Module", envir = .GlobalEnv)) {
  Rcpp::loadModule("DUSTMODULEMD", TRUE)
}

#' Run the MD Dust Change Point Detection
#'
#' This function performs change point detection on multi-dimensional data using the DUST algorithm.
#'
#' @param data A numeric matrix. The time series data on which change point detection is performed. Each row represents a time-series.
#' @param penalty A numeric value. The penalty applied for adding a new change point. By default, it is set to \code{2 x nrow(data) x log(length(data))}.
#' @param constraints_l An integer. The number of left constraints to be considered in the DUST pruning test.
#' @param constraints_r An integer. The number of right constraints to be considered in the DUST pruning test.
#' @param model A character string. Specifies the model used for change point detection. Default is \code{"gauss"}. Possible values could include \code{"gauss"}, \code{"poisson"}, \code{"exp"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}, \code{"variance"}.
#' @param method A character string specifying the method used to handle indices and pruning tests in the algorithm. The default is \code{"detIndex_Eval4"}, which is the fastest method for the chosen model. Other available methods are:
#' \itemize{
#'   \item \code{"randIndex_Eval0"} to \code{"randIndex_Eval6"}: Random index-based methods with different dual maximization algorithm (0 through 6).
#'   \item \code{"detIndex_Eval0"} to \code{"detIndex_Eval6"}: Deterministic index-based methods  with different dual maximization algorithm (0 through 6).
#' }
#' Here are the current available algorithms (\code{Eval4} is often the most efficient one)
#' \itemize{
#'   \item \code{"Eval0"}: random evaluation of the dual (with uniform distribution)
#'   \item \code{"Eval1"}:
#'   \item \code{"Eval2"}:
#'   \item \code{"Eval3"}:
#'   \item \code{"Eval4"}:
#'   \item \code{"Eval5"}:
#'   \item \code{"Eval6"}:
#' }
#' @param nbLoops An integer. The number of loops to run in the max dual optimization algorithm. Default is 10.
#'
#' @return A list containing the information computed by the DUST algorithm.
#' \itemize{
#'   \item \code{changepoints}: the sequence of optimal change points solving our penalized optimization problem
#'   \item \code{lastIndexSet}: the last non-pruned indices at time step n (= data length)
#'   \item \code{nb}: vector or size n (= data length) recording the number of non-pruned indices over time
#'   \item \code{costQ}: vector or size n (= data length) recording the optimal (penalized) segmentation cost over time
#' }
#'
#' @examples
#' mu <- c(rep(0,100), rep(1, 100))
#' data <- matrix(rnorm(100 * 2, mu), nrow = 2, ncol = 100, byrow = FALSE)
#' # dust.MD(data)
#' # param <- data.frame(c(1,1,1), c(2,1,1),c(1,0.1,2))
#' # data <- dataGenerator_MD(chpts = c(40,60,100), parameters = param, type = "exp")
#' # dust.MD(data, penalty =  6*log(100), model = "exp", method = "detIndex_Eval4")
#' @export
dust.MD <- function(
    data = data
    , penalty = 2*nrow(data)*log(length(data))
    , constraints_l = nrow(data)
    , constraints_r = 0
    , model = "gauss"
    , method = "detIndex_Eval4"
    , nbLoops = 10
)
{
  if (!is.null(nbLoops) && nbLoops < 0)
  {
    stop("negative nbLoops passed to dust.MD")
  }
  partitioner <- new(DUST_MD, model, method, nbLoops)
  return(partitioner$one_dust(data, penalty, constraints_l, constraints_r))
}




