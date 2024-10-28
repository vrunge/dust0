source("Simulations/Section_Simulations/ENGINE.R")

SIMULATION.PARAMETERS = list(
  MODELS = model.compare,
  METHOD = "", # "" -> fastest
  SIZES = c(1e4, 1e8),
  N.CORES = 20, # number of cores used simultaneously (at most)
  N.JOBS = 100, # number of executions per configuration
  FILE.NAME = "nb8" # simulation will be saved as FILE.NAME_CURRENT.MODEL.csv
)


# ----------------------- #
# --- ############### --- #
# --- ### NB PLOT ### --- #
# --- ############### --- #
# ----------------------- #

nb.plot = function(i, engine, sizes, data.generator, beta.generator)
  do.call(
    rbind,
    lapply(
      sizes,
      function(size)
      {
        engine@size = size
        engine@beta = beta.generator(size)

        if (size > 1e4) { ftr = rep(c(FALSE, TRUE), c(size %/% 1e4 - 1, 1)) }
        else { ftr = TRUE }

        d = get_param(engine)
        d[['exec']] = i
        d[['time']] = (1:size)[ftr]
        ex = compute(engine, data.generator(size))
        d[['nb.cpts']] = length(ex$changepoints)
        d[['nb']] = ex$nb[ftr]
        as.data.frame.list(d)
      }
    )
  )


# --- ### OUTPUT ### --- #

beta.generator = function(n, d = 1) 2 * d * log(n)

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
