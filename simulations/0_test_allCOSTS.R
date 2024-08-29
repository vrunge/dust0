


n <- 10^3
beta <- 2 * log(n)

type = "gauss"
param = 0.5
myData <- dataGenerator_1D(n, parameters = param, type = type)
d <- dust.partitioner.1D(model = type, method = "randIndex_randEval", nbLoops = 1)

d$quick(data = myData, penalty = beta)



####

n =3*10^5
beta = 2*log(n)
y <- dataGenerator_1D(chpts = n, parameters = 1, type = "poisson")
t1 <- system.time(dust.partitioner.1D(model = "poisson", method = "fastest")$quick(data = y, penalty = beta))
t2 <- system.time(dust.partitioner.1D(model = "poisson", method = "randIndex_randEval")$quick(data = y, penalty = beta))
t1
t2


n = 10^7
beta = 2*log(n)
y <- dataGenerator_1D(chpts = n, parameters = 1, type = "poisson")
t2 <- system.time(dust.partitioner.1D(model = "poisson", method = "fastest")$quick(data = y, penalty = beta))
t2







library(fpopw)
n = 10^7
beta = 2*log(n)/50
y <- dataGenerator_1D(chpts = n, parameters = 0, type = "gauss")
t3 <- system.time(dust.partitioner.1D(model = "gauss", method = "fastest")$quick(data = y, penalty = beta))
t3
s <- system.time(fpopw::Fpop(y, beta))
s
(t3[[1]] - s[[1]])/s[[1]]


n = 10^8
y <- dataGenerator_1D(chpts = n, parameters = 0, type = "gauss")
t3 <- system.time(dust.partitioner.1D(model = "gauss", method = "fastest")$quick(data = y, penalty = beta))
t3
s <- system.time(fpopw::Fpop(y, beta))
s
(t3[[1]] - s[[1]])/s[[1]]


n = 10^7
beta = 0
y <- dataGenerator_1D(chpts = n, parameters = 0, type = "gauss")
res3 <- dust.partitioner.1D(model = "gauss", method = "fastest")$quick(data = y, penalty = beta)
length(res3$changepoints)

v <- dust_R_1D(y,penalty = beta,
               type = "gauss",
               pruningOpt = 0)

res3$changepoints
v$changepoints
res3$costQ - v$costQ


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


##################################################################
##################################################################


nb <- 9
n <- 10^6
beta <- 2 * log(n)

test <- function(n = 10^4, type = "gauss", param = 0.5, beta = 0)
{
  if(beta == 0){beta <- 2 * log(n)}
  myData <- dataGenerator_1D(n, parameters = param, type = type)

  if(type == "binom"){myData <- myData/max(myData)}
  t1 <- system.time(dust.partitioner.1D(model = type, method = "randIndex_randEval", nbLoops = nb)$quick(data = myData, penalty = beta))
  t2 <- system.time(dust.partitioner.1D(model = type, method = "randIndex_detEval", nbLoops = nb)$quick(data = myData, penalty = beta))
  t3 <- system.time(dust.partitioner.1D(model = type, method = "fastest", nbLoops = nb)$quick(data = myData, penalty = beta))
  u <- dust.partitioner.1D(model = type, method = "fastest")$quick(data = myData, penalty = beta)
  print(u$changepoints)
  #v <- dust_R_1D(myData, beta, type = type, pruningOpt = 2)
  #print(u$changepoints)
  #print(v$changepoints)
  return(c(t1[[1]],t2[[1]],t3[[1]]))
}

test(n, "gauss", beta = beta)
test(n, "poisson", beta = beta)
test(n, "exp", beta = beta)
test(n, "geom", beta = beta)
test(n, "bern", beta = beta)
test(n, "binom", beta = beta)
test(n, "negbin", beta = beta)

##################################################################
##################################################################



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
