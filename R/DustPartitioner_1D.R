
## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##

Rcpp::loadModule("DUSTMODULE1D", TRUE)

Rcpp::loadModule("FLATDUST1D", TRUE) #with reserve(n + 1);
Rcpp::loadModule("FLAT2DUST1D", TRUE) #with [n+1]

Rcpp::loadModule("FLATOP1D", TRUE)

## --------------------------------- ##
## ----///////////////////////// --- ##
## --- //    1D DUST object   // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.object.1D
#'
#' @description Generates a DUST 1D object
#'
#' @param model A character string specifying the model for the data. The default is \code{gauss}. Available models are:
#' \itemize{
#'   \item \code{"gauss"}: Assumes the data follows a Gaussian distribution with known variance
#'   \item \code{"poisson"}: Assumes the data follows a Poisson distribution, typically for count data.
#'   \item \code{"exp"}: Assumes the data follows an exponential distribution.
#'   \item \code{"geom"}: Assumes the data follows a geometric distribution.
#'   \item \code{"bern"}: Assumes the data follows a Bernoulli distribution, typically for binary data.
#'   \item \code{"binom"}: Assumes the data follows a Binomial distribution, for experiments with a fixed number of trials.
#'   \item \code{"negbin"}: Assumes the data follows a Negative Binomial distribution, for overdispersed count data.
#'   \item \code{"variance"}: Assumes the data follows a Gaussian distribution with unknown variance and null mean.
#' }
#' @param method A character string specifying the method used to handle indices and pruning tests in the algorithm. The default is \code{detIndex_Eval4}, which automatically selects the fastest method for the chosen model. Other available methods are:
#' \itemize{
#'   \item \code{"randIndex_Eval0"} to \code{"randIndex_Eval6"}: Random index-based methods with different dual maximization algorithm (0 through 5).
#'   \item \code{"detIndex_Eval0"} to \code{"detIndex_Eval6"}: Deterministic index-based methods  with different dual maximization algorithm (0 through 5).
#' }
#' Here are the current available algorithms (\code{Eval4} is often the most efficient one)
#' \itemize{
#'   \item \code{"Eval0"}: random evaluation of the dual (with uniform distribution)
#'   \item \code{"Eval1"}: max value with closed formula (gauss model only), otherwise no pruning performed and we get the (slow) OP algorithm
#'   \item \code{"Eval2"}: golden-section search.
#'   \item \code{"Eval3"}: binary search. At each step, we evaluate the tangent line to the current point at its max to stop the search at early step (when possible)
#'   \item \code{"Eval4"}: quasi-Newton method with armijo condition
#'   \item \code{"Eval5"}: PELT rule
#'   \item \code{"Eval6"}: OP rule
#' }
#' @param nbLoops number of iterations in the algorithm for maximizing the dual function
#'
#' @return a DUST 1D object that provides methods:
#' \itemize{
#'   \item \code{append}, for adding new data to the vector of data to be analysed;
#'   \item \code{update_partition}, for update the optimal partition with new added data using append method;
#'   \item \code{get_partition}, for retrieving the optimal partition once it has been computed.
#' }
#' @examples
#' dust.object.1D()
dust.object.1D <- function(
    model = "gauss"
    , method = "detIndex_Eval4"
    , nbLoops = 10
)
{
  object <- new(DUST_1D, model, method, nbLoops)
  assign(
    "append",
    function(data, penalty = NULL)
      object$append_c(data, penalty),
    envir = object
  )
  return(object)
}

