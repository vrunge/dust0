
## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##

Rcpp::loadModule("DUSTMODULE1D", TRUE)


## --------------------------------- ##
## ----///////////////////////// --- ##
## --- // 1D DUST Partitioner // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.partitioner.1D
#'
#' @description Generates a DUST partitioner
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
#' }
#' @param alpha controls the randomness of the random methods. for computational efficiency purposes, a vector of random values is generated upon initializing the partitioner object.
#' the lower the value of alpha, the larger this vector, meaning a selection of indices closer to "true" randomness.
#' @param nbLoops number of iterations in the algorithm for maximizing the dual function
#'
#' @return a DUST partitioner object that provides methods:
#' \itemize{
#'   \item \\code{fit}, for fitting the data;
#'   \item \code{compute}, once fit has been called, for computing the optimal partition of the data;
#'   \item \code{get_partition}, for retrieving the optimal partition once it has been computed; and
#'   \item \code{quick}, a wrapper for the 3 methods.
#' }
#' @examples
#' dust.partitioner.1D()
dust.partitioner.1D <- function(
    model = "gauss"
    , method = "fastest"
    , alpha = 1e-9
    , nbLoops = 10
)
{
  partitioner <- new(DUST_1D, model, method, alpha, nbLoops)

  assign(
    "init",
    function(data, penalty = NULL)
      partitioner$init_raw(data, penalty),
    envir = partitioner
  )

  assign(
    "quick",
    function(data, penalty = NULL)
      partitioner$quick_raw(data, penalty),
    envir = partitioner
  )

  return(partitioner)
}



####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##

Rcpp::loadModule("DUSTMODULEMeanVar", TRUE)


## --------------------------------- ##
## ----///////////////////////// --- ##
## --- // 2D DUST Partitioner // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.partitioner.meanVar
#'
#' @description Generates a DUST partitioner
#' @param method the method for handling the indices and pruning tests in the algorithm. defaults to "fastest", which returns the fastest method available on input model. available methods: "randIndex_randEval" prunes indices based on a random dual test, where the value of the dual is evaluated at a random point, and the dual is defined by a random index; "randIndex_detEval", where the dual is defined by a random index but evaluated at its maximum.
#' @param alpha controls the randomness of the random methods. for computational efficiency purposes, a vector of random values is generated upon initializing the partitioner object. the lower the value of alpha, the larger this vector, meaning a selection of indices closer to "true" randomness.
#' @param nbLoops number of iteration in the optimization algorithm for maximizing the dual function
#' @return a DUST partitioner object that provides methods : fit, for fitting the data; compute, once fit has been called, for computing the optimal partition of the data; get_partition, for retrieving the optimal partition once it has been computed; and quick, a wrapper for the other 3 methods.
#' @examples
#' dust.partitioner.meanVar()
dust.partitioner.meanVar <- function(
     method = "fastest"
    , alpha = 1e-9
    , nbLoops = 10
)
{
  partitioner <- new(DUST_meanVar, method, alpha, nbLoops)

  assign(
    "init",
    function(data, penalty = NULL)
      partitioner$init_raw(data, penalty),
    envir = partitioner
  )

  assign(
    "quick",
    function(data, penalty = NULL)
      partitioner$quick_raw(data, penalty),
    envir = partitioner
  )

  return(partitioner)
}


####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################

## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##

Rcpp::loadModule("DUSTMODULEreg", TRUE)


## --------------------------------- ##
## ----///////////////////////// --- ##
## --- // 2D DUST Partitioner // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.partitioner.reg
#'
#' @description Generates a DUST partitioner
#' @param method the method for handling the indices and pruning tests in the algorithm.
#' defaults to "fastest", which returns the fastest method available on input model.
#' available methods: "randIndex_randEval" prunes indices based on a random dual test,
#' where the value of the dual is evaluated at a random point, and the dual is defined
#' by a random index; "randIndex_detEval", where the dual is defined by a random index
#' but evaluated at its maximum.
#' @param alpha controls the randomness of the random methods. for computational efficiency purposes, a vector of random values is generated upon initializing the partitioner object. the lower the value of alpha, the larger this vector, meaning a selection of indices closer to "true" randomness.
#' @param nbLoops number of iteration in the optimization algorithm for maximizing the dual function
#' @return a DUST partitioner object that provides methods : fit, for fitting the data; compute, once fit has been called, for computing the optimal partition of the data; get_partition, for retrieving the optimal partition once it has been computed; and quick, a wrapper for the other 3 methods.
#' @examples
#' dust.partitioner.reg()
dust.partitioner.reg <- function(
    method = "fastest"
    , alpha = 1e-9
    , nbLoops = 10
)
{
  partitioner <- new(DUST_reg, method, alpha, nbLoops)

  assign(
    "init",
    function(data, penalty = NULL)
      partitioner$init_raw(data, penalty),
    envir = partitioner
  )

  assign(
    "quick",
    function(data, penalty = NULL)
      partitioner$quick_raw(data, penalty),
    envir = partitioner
  )

  return(partitioner)
}


## ----------------------------------- ##
## --- /////////////////////////// --- ##
## --- // Importing C++ Modules // --- ##
## --- /////////////////////////// --- ##
## ----------------------------------- ##

Rcpp::loadModule("DUSTMODULEMD", TRUE)


## --------------------------------- ##
## ----///////////////////////// --- ##
## --- // MD DUST Partitioner // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.partitioner.MD
#'
#' @description Generates a DUST partitioner
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
#' }
#' @param alpha controls the randomness of the random methods. for computational efficiency purposes, a vector of random values is generated upon initializing the partitioner object.
#' the lower the value of alpha, the larger this vector, meaning a selection of indices closer to "true" randomness.
#' @param nbLoops number of iterations in the algorithm for maximizing the dual function
#'
#' @return a DUST partitioner object that provides methods:
#' \itemize{
#'   \item \\code{fit}, for fitting the data;
#'   \item \code{compute}, once fit has been called, for computing the optimal partition of the data;
#'   \item \code{get_partition}, for retrieving the optimal partition once it has been computed; and
#'   \item \code{quick}, a wrapper for the 3 methods.
#' }
#' @examples
#' dust.partitioner.1D()
dust.partitioner.MD <- function(
    model = "gauss"
    , method = "fastest"
    , alpha = 1e-9
    , nbLoops = 10
)
{
  partitioner <- new(DUST_MD, model, method, alpha, nbLoops)

  assign(
    "init",
    function(data, penalty = NULL, nb_max = NULL)
      partitioner$init_raw(data, penalty, nb_max),
    envir = partitioner
  )

  assign(
    "quick",
    function(data, penalty = NULL, nb_max = NULL)
      partitioner$quick_raw(data, penalty),
    envir = partitioner
  )

  return(partitioner)
}

