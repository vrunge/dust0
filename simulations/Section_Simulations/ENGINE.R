# --- ### GENERATORS ### --- #

# devtools::install_github("vrunge/dust")

library(dust0)

# Model names
model.names = c("gauss", "poisson", "exp", "geom", "bern", "negbin", "variance")
model.compare = c("gauss", "poisson", "exp", "negbin", "variance")

base.parameters = c(
  "gauss" = 0,
  "poisson" = 3,
  "exp" = 1,
  "geom" = .5,
  "bern" = .5,
  "negbin" = .5,
  "variance" = 1
)

# Data generators, adjusted for a certain behaviour under 2 * log(n) penalty
data.generators = list(
  "gauss" = function(size, parameters) rnorm(sum(size), rep(parameters, size)),
  "poisson" = function(size, parameters) rpois(sum(size), rep(parameters, size)),
  "exp" = function(size, parameters) rexp(sum(size), rep(parameters, size)),
  "geom" = function(size, parameters) dataGenerator_1D(cumsum(size), parameters = parameters, type = "geom"),
  "bern" = function(size, parameters) dataGenerator_1D(cumsum(size), parameters = parameters, type = "bern"),
  "negbin" = function(size, parameters) dataGenerator_1D(cumsum(size), parameters = parameters, type = "negbin"),
  "variance" = function(size, parameters) dataGenerator_meanVar(cumsum(size), means = rep(0, length(parameters)), sds = parameters)
)

f.generators = function(model) function(size, parameters = base.parameters[model]) data_normalization_1D(data.generators[[model]] (size, parameters), model)

model.mult = list(
  gauss = 1,
  exp = 1,
  poisson = 1/3,
  negbin = .1,
  geom = 1,
  bern = 1,
  variance = 1
)

data.parameters = list(
  "gauss" = c(0, 1),
  "poisson" = c(3, 5),
  "exp" = c(1, 1/2),
  "geom" = c(.5, .7),
  "bern" = c(.5, .7),
  "negbin" = c(.5, .7),
  "variance" = c(1, 2)
)

cpt.generators = function(density = 0) function(model = "gauss")
{
  generator = data.generators[[model]]
  lapply(
    density,
    function(d)
      list(
        density = d,
        generator = function(size)
        {
          if (d > 0)
          {
            n = c(rep(d, size %/% d))
            if (size %% d > 0)
              n = c(n, size %% d)
          } else n = size
          parameters = rep(data.parameters[[model]], length.out = length(n))
          return(generator(n, parameters))
        }
      )
  )
}

# ----------------------------- #
# --- ##################### --- #
# --- # Simulation engine # --- #
# --- ##################### --- #
# ----------------------------- #

setClass(

  "simulation.engine",

  slots = c(
    algo = "character"
    , model = "character"
    , method = "character"
    , size = "numeric"
    , cpt.density = "numeric"
    , beta = "numeric"
    , alpha = "numeric"
    , nbLoops = "numeric"
    , m1 = "numeric"
  )

  , prototype = list(
    algo = "dust"
    , model = "gauss"
    , method = ""
    , size = 0
    , cpt.density = 0
    , beta = NA_real_
    , alpha = 1e-9
    , nbLoops = 10
    , m1 = 1e-2
  )

)

setGeneric("compute", function(x, data) standardGeneric("compute"))
setGeneric("get_param", function(x) standardGeneric("get_param"))
setGeneric("copy_engine", function(x) standardGeneric("copy_engine"))

setMethod(

  "compute"

  , signature = c(x = "simulation.engine", data = "numeric")

  , function(x, data)
  {

    beta = x@beta
    if (x@algo == "dust") {

      model = x@model

      if (x@method %in% c("", "fastest")) {

        method = ifelse(model == "gauss", "detIndex_Eval1", "detIndex_Eval4")

      }

      else method = method

      beta = beta * model.mult[[model]]

      return(dust.1D(data, beta, model = model, method = method, x@nbLoops, x@alpha))

    } else if (x@algo == "fpop") {

      if (x@model %in% c("", "gauss")) {

        return(fpopw::Fpop(data, beta))

      } else if (!(x@model %in% model.compare)) {

        stop(x@model, " model invalid for FPOP comparison")

      } else {

        model = x@model

        beta = beta * model.mult[[model]]

        myGraph = gfpop::graph(type = "std", penalty = 2 * beta)

        return(gfpop::gfpop(data, myGraph, type = model))

      }

    }

  }

)

setMethod(

  "get_param"

  , signature = c(x = "simulation.engine")

  , function(x)
  {

    model = ifelse(x@model == "", "gauss", x@model)
    method = x@method

    if (x@algo == "dust") {

      if (method %in% c("", "fastest"))
        method = ifelse(model == "gauss", "detIndex_Eval1", "detIndex_Eval4")

    }

    else method = ""

    list(
      algo = x@algo
      , model = model
      , method = method
      , size = x@size
      , cpt.density = x@cpt.density
      , beta = x@beta * model.mult[[model]]
      , alpha = x@alpha
      , nbLoops = x@nbLoops
      , m1 = x@m1
    )

  }

)

setMethod(

  "copy_engine"

  , signature = c(x = "simulation.engine")

  , function(x)
  {
    out = new("simulation.engine")
    out@algo = x@algo
    out@model = x@model
    out@method = x@method
    out@size = x@size
    out@cpt.density = x@cpt.density
    out@beta = x@beta
    out@alpha = x@alpha
    out@nbLoops = x@nbLoops
    out@nbLoops = x@nbLoops
    return(out)
  }

)


# --- ### FUNCTION ### --- #

simulate = function(engine, job.generator, models, method, sizes, beta.generator, data.generators, n.cores, n.jobs, fname = "simu")
{
  print(paste("Executions per setup:", n.jobs))

  do.call(
    rbind
    , lapply(
      models,
      function(name)
      {

        print(paste("generating simulations for", name, "model"))

        engine@model = name
        generator = data.generators(name)
        out = do.call(
          rbind,
          parallel::mclapply(
            # lapply(
            1:n.jobs
            , job.generator
            , engine = copy_engine(engine)
            , sizes = sizes
            , data.generator = generator
            , beta.generator = beta.generator
            , mc.cores = n.cores
            , mc.preschedule = FALSE
          )
        )
        path <- "Simulations/Section_Simulations/"
        write.csv(out, paste0(path, paste(fname, name, sep = "_"), '.csv'))

        return(out)
      }
    )
  )
}
