
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
#' @param model the underlying model of the data. defaults to "gauss". available models: c("gauss", "poisson)
#' @param method the method for handling the indices and pruning tests in the algorithm. defaults to "fastest", which returns the fastest method available on input model. available methods: "randIndex_randEval" prunes indices based on a random dual test, where the value of the dual is evaluated at a random point, and the dual is defined by a random index; "randIndex_detEval", where the dual is defined by a random index but evaluated at its maximum.
#' @param alpha controls the randomness of the random methods. for computational efficiency purposes, a vector of random values is generated upon initializing the partitioner object. the lower the value of alpha, the larger this vector, meaning a selection of indices closer to "true" randomness.
#' @return a DUST partitioner object that provides methods : fit, for fitting the data; compute, once fit has been called, for computing the optimal partition of the data; get_partition, for retrieving the optimal partition once it has been computed; and quick, a wrapper for the other 3 methods.
#' @examples

dust.partitioner.1D = function(
    model = "gauss"
    , method = "fastest"
    , alpha = 1e-9
)
{
  partitioner = new(DUST_1D, model, method, alpha)

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

