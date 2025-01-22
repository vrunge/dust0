
Rcpp::loadModule("DUSTMODULEMeanVar", TRUE)

## --------------------------------- ##
## ----///////////////////////// --- ##
## --- // 2D DUST Partitioner // --- ##
## ----///////////////////////// --- ##
## --------------------------------- ##

#' dust.partitioner.meanVar
#'
#' @description Generates a DUST partitioner
#' @param method the method for handling the indices and pruning tests in the algorithm. defaults to "detIndex_Eval4", which returns the fastest method available on input model. available methods: "randIndex_randEval" prunes indices based on a random dual test, where the value of the dual is evaluated at a random point, and the dual is defined by a random index; "randIndex_detEval", where the dual is defined by a random index but evaluated at its maximum.
#' @param nbLoops number of iteration in the optimization algorithm for maximizing the dual function
#' @return a DUST partitioner object that provides methods : fit, for fitting the data; compute, once fit has been called, for computing the optimal partition of the data; get_partition, for retrieving the optimal partition once it has been computed; and quick, a wrapper for the other 3 methods.
#' @examples
#' dust.partitioner.meanVar()
dust.partitioner.meanVar <- function(
    method = "detIndex_Eval4"
    , nbLoops = 10
)
{
  partitioner <- new(DUST_meanVar, method, nbLoops)

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
#' defaults to "detIndex_Eval4", which returns the fastest method available on input model.
#' available methods: "randIndex_randEval" prunes indices based on a random dual test,
#' where the value of the dual is evaluated at a random point, and the dual is defined
#' by a random index; "randIndex_detEval", where the dual is defined by a random index
#' but evaluated at its maximum.
#' @param nbLoops number of iteration in the optimization algorithm for maximizing the dual function
#' @return a DUST partitioner object that provides methods : fit, for fitting the data; compute, once fit has been called, for computing the optimal partition of the data; get_partition, for retrieving the optimal partition once it has been computed; and quick, a wrapper for the other 3 methods.
#' @examples
#' dust.partitioner.reg()
dust.partitioner.reg <- function(
    method = "detIndex_Eval4"
    , nbLoops = 10
)
{
  partitioner <- new(DUST_reg, method, nbLoops)

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
