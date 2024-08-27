

n <- 10^4
beta <- 5 * log(n)

test <- function(n = 10^4, type = "gauss", param = 0.5)
{
  myData <- dataGenerator_1D(n, parameters = param, type = type)

  t1 <- system.time(dust.partitioner.1D(model = type, method = "randIndex_randEval")$quick(data = myData, penalty = beta))
  t2 <- system.time(dust.partitioner.1D(model = type, method = "randIndex_detEval")$quick(data = myData, penalty = beta))
  t3 <- system.time(dust.partitioner.1D(model = type, method = "fastest")$quick(data = myData, penalty = beta))
  u <- dust.partitioner.1D(model = type, method = "randIndex_randEval")$quick(data = myData, penalty = beta)
  v <- dust_R_1D(myData, beta, type = type)
  print(u$changepoints)
  print(v$changepoints)
  return(c(t1[[1]],t2[[1]],t3[[1]]))
}


test(n, "gauss")
test(n, "poisson")
test(n, "exp")
test(n, "geom")
test(n, "bern")

type <- "bern"
param <- 0.5
myData <- dataGenerator_1D(n, parameters = param, type = type)

d <- dust.partitioner.1D(model = type, method = "randIndex_randEval")$quick(data = myData, penalty = beta)
d$changepoints
dust.partitioner.1D(model = type, method = "randIndex_detEval")$quick(data = myData, penalty = beta)
dust.partitioner.1D(model = type, method = "fastest")$quick(data = myData, penalty = beta)



test(n, "binom")
test(n, "negbin")

################################################################################
################################################################################


n <- 1000
beta <- 2 * log(n)
chgtpt <- c(0.1, 0.3, 0.5, 0.8, 1)
myData <- dataGenerator(n, chgtpt, c(10, 1, 10, 1, 10), type = "poisson")
#plot(myData)

GRAPH <- graph(type = "std", penalty = beta)

resGFPOP <- gfpop(data = myData, type = "poisson", mygraph = GRAPH)
resDUST <- dust.partitioner.1D(model = "poisson", method = "randIndex_randEval")$quick(data = myData, penalty = beta)


resGFPOP$globalCost
