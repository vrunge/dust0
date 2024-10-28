source("scripts/ENGINE.R")

SIMULATION.PARAMETERS = list(
  MODELS = model.compare,
  METHOD = "", # "" -> fastest
  SIZES = c(1e3, 1e6),
  N.CORES = 20, # number of cores used simultaneously (at most)
  N.JOBS = 100, # number of executions per configuration
  FILE.NAME = "beta6" # simulation will be saved as FILE.NAME_CURRENT.MODEL.csv
)

SIMULATION.PARAMETERS[["BETA.VALUES"]] = exp(seq(log(1e-3), log(20), length.out = 21))


# ------------------------- #
# --- ################# --- #
# --- ### BETA TIME ### --- #
# --- ################# --- #
# ------------------------- #

beta.time = function(i, engine, sizes, data.generator, beta.generator)
  do.call(
    rbind,
    lapply(
      sizes,
      function(size)
      {
        beta.v = beta.generator(size)
        
        do.call(
          rbind,
          lapply(
            beta.v,
            function(beta)
            {
              
              engine@beta = beta
              
              y = data.generator(size)
              
              dust = get_param(engine)
              dust[['exec']] = i
              dust[['size']] = size
              dust[['exec.time']] = microbenchmark::microbenchmark(ex <- compute(engine, y), times = 1)$time[1]
              dust[['nb.cpts']] = length(ex$changepoints)
              dust[['nb']] = tail(ex$nb, 1)
              
              rm(ex)
              
              as.data.frame.list(dust)
              
            }
          )
        )
        
      }
    )
  )


# --- ### OUTPUT ### --- #

beta.generator = function(n, d = 1) SIMULATION.PARAMETERS$BETA.VALUES * d * log(n)

out <- simulate(
  engine = new("simulation.engine", algo = "dust")
  , job.generator = nb.plot
  , models = SIMULATION.PARAMETERS$MODELS
  , method = SIMULATION.PARAMETERS$METHOD
  , sizes = SIMULATION.PARAMETERS$SIZES
  , beta.generator = beta.generator
  , data.generators = f.generators
  , n.cores = SIMULATION.PARAMETERS$N.CORES
  , n.jobs = SIMULATION.PARAMETERS$N.JOBS
  , fname = SIMULATION.PARAMETERS$FILE.NAME
)
