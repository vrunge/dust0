source("scripts/ENGINE.R")

SIMULATION.PARAMETERS = list(
  MODELS = model.compare,
  METHOD = "", # "" -> fastest
  SIZES = c(1e3, 1e6),
  N.CORES = 20, # number of cores used simultaneously (at most)
  N.JOBS = 10, # number of executions per configuration
  FILE.NAME = "cpt6" # simulation will be saved as FILE.NAME_CURRENT.MODEL.csv
)

SIMULATION.PARAMETERS[["DENSITY.VALUES"]] = seq(0, 50, length.out = 6)


# ------------------------ #
# --- ################ --- #
# --- ### CPT TIME ### --- #
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

        do.call(

          rbind,

          lapply(

            data.generator,

            function(density)
            {

              engine@cpt.density = density$density

              y = density$generator(size)

              engine@algo = "dust"
              dust = get_param(engine)
              dust[['exec']] = i
              dust[['exec.time']] = microbenchmark::microbenchmark(ex <- compute(engine, y), times = 1)$time[1]
              dust[['nb.cpts']] = length(ex$changepoints)
              dust[['nb']] = tail(ex$nb, 1)

              engine@algo = "fpop"
              fpop = get_param(engine)
              fpop[['exec']] = i
              fpop[['exec.time']] = microbenchmark::microbenchmark(ex <- compute(engine, y), times = 1)$time[1]
              fpop[['nb.cpts']] = length(ex$changepoints) + (engine@model == "gauss")
              fpop[['nb']] = NA

              rbind(dust, fpop)

            }

          )

        )

      }

    )

  )

# --- ### OUTPUT ### --- #

beta.generator = function(n, d = 1) 2 * d * log(n)

out <- simulate(
  engine = new("simulation.engine")
  , job.generator = log.time.compare
  , models = model.compare
  , method = ""
  , sizes = exp(seq(log(1e2), log(5e2), length.out = 3))
  , beta.generator = beta.generator
  , data.generators = cpt.generators(SIMULATION.PARAMETERS$DENSITY.VALUES)
  , n.cores = SIMULATION.PARAMETERS$N.CORES
  , n.jobs = SIMULATION.PARAMETERS$N.JOBS
  , fname = SIMULATION.PARAMETERS$FILE.NAME
)
