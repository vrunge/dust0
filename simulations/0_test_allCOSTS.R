

n = 10^7
beta = 2*log(n)
y <- dataGenerator_1D(chpts = n, parameters = 0, type = "gauss")
t3 <- system.time(dust.partitioner.1D(model = "gauss", method = "fastest")$quick(data = y, penalty = beta))
t3

n = 10^8
beta = 2*log(n)
y <- dataGenerator_1D(chpts = n, parameters = 0, type = "gauss")
t3 <- system.time(dust.partitioner.1D(model = "gauss", method = "fastest")$quick(data = y, penalty = beta))
t3




#res1 <- dust.partitioner.1D(model = "gauss", method = "randIndex_randEval")$quick(data = y, penalty = beta)
res3 <- dust.partitioner.1D(model = "gauss", method = "fastest")$quick(data = y, penalty = beta)

res1$changepoints
res3$changepoints


###################################################################################################


y
res2$costQ

t1$changepoints
t2$changepoints
t3$changepoints



n <- 10^2
beta <- 2 * log(n)/1000
type = "bern"
param = c(0.2,0.6)

myData <- dataGenerator_1D(c(n/2,n), parameters = param, type = type)
u <- dust.partitioner.1D(model = type, method = "fastest")$quick(data = myData, penalty = beta)
v <- dust_R_1D(myData,
               beta,
               type = type,
               pruningOpt = 0)

print(u$changepoints)
print(v$changepoints)
v


n <- 10^4
beta <- 5 * log(n)/10

test <- function(n = 10^4, type = "gauss", param = 0.5)
{
  myData <- dataGenerator_1D(n, parameters = param, type = type)

  if(type == "binom"){myData <- myData/max(myData)}
  t1 <- system.time(dust.partitioner.1D(model = type, method = "randIndex_randEval")$quick(data = myData, penalty = beta))
  t2 <- system.time(dust.partitioner.1D(model = type, method = "randIndex_detEval")$quick(data = myData, penalty = beta))
  t3 <- system.time(dust.partitioner.1D(model = type, method = "fastest")$quick(data = myData, penalty = beta))
  u <- dust.partitioner.1D(model = type, method = "randIndex_randEval")$quick(data = myData, penalty = beta)
  v <- dust_R_1D(myData, beta, type = type, pruningOpt = 2)
  print(u$changepoints)
  print(v$changepoints)
  return(c(t1[[1]],t2[[1]],t3[[1]]))
}


test(n, "gauss")
test(n, "poisson")
test(n, "exp")
test(n, "geom")
test(n, "bern") ### PB
test(n, "binom")
test(n, "negbin")

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
