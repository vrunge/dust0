source("scripts/ENGINE.R")

SIMULATION.PARAMETERS = list(
  MODELS = model.compare,
  METHOD = "", # "" -> fastest
  SIZES = exp(seq(log(100), log(1e6), length.out = 100)),
  N.CORES = 20, # number of cores used simultaneously (at most)
  N.JOBS = 100, # number of executions per configuration
  FILE.NAME = "reg62"
)

# ------------------------ #
# --- ################ --- #
# --- ### LOG TIME ### --- #
# --- ################ --- #
# ------------------------ #

log.time.compare = function(i, engine, sizes, data.generator, beta.generator)
  do.call(
    rbind,
    lapply(
      sizes,
      function(size)
      {
        engine@size = size
        engine@beta = beta.generator(size)

        y = data.generator(size)

        engine@algo = "dust"
        dust = get_param(engine)
        dust[['exec']] = i
        dust[['size']] = size
        dust[['exec.time']] = microbenchmark::microbenchmark(ex <- compute(engine, y), times = 1)$time[1]
        dust[['nb.cpts']] = length(ex$changepoints)
        dust[['nb']] = tail(ex$nb, 1)

        engine@algo = "fpop"
        fpop = get_param(engine)
        fpop[['exec']] = i
        fpop[['size']] = size
        fpop[['exec.time']] = microbenchmark::microbenchmark(ex <- compute(engine, y), times = 1)$time[1]
        fpop[['nb.cpts']] = length(ex$changepoints) + (engine@model == "gauss")
        fpop[['nb']] = NA

        rm(ex)

        rbind(dust, fpop)

      }
    )
  )


# --- ### OUTPUT ### --- #

beta.generator = function(n, d = 1) 2 * d * log(n)

out <- simulate(
  engine = new("simulation.engine")
  , job.generator = log.time.compare
  , models = SIMULATION.PARAMETERS$MODELS
  , method = SIMULATION.PARAMETERS$METHOD
  , sizes = SIMULATION.PARAMETERS$SIZES
  , beta.generator = beta.generator
  , data.generators = f.generators
  , n.cores = SIMULATION.PARAMETERS$N.CORES
  , n.jobs = SIMULATION.PARAMETERS$N.JOBS
  , fname = SIMULATION.PARAMETERS$FILE.NAME
)
