library(dust0)


tpe <- "poisson"
n <- 10^2
beta <- 2*log(n)
y <- dataGenerator_1D(n, parameters = 1, type = tpe)
plot_dual_1D(data = y, s1 = 95, s2 = 90, type = tpe, mu = 1:9999/10000, penalty = 3, OP = T)




te <- "randIndex_Eval0"

n <- 10^3
beta <- 2*log(n)
y <- rnorm(n)
res1 <- dust.partitioner.1D(method = "detIndex_Eval0")$quick(data = y, penalty = beta)
res2 <- dust.partitioner.1D(method = "detIndex_Eval1")$quick(data = y, penalty = beta)
res3 <- dust.partitioner.1D(method = "detIndex_Eval2")$quick(data = y, penalty = beta)
res4 <- dust.partitioner.1D(method = "detIndex_Eval3")$quick(data = y, penalty = beta)

plot(res1$nb, type = 'l')
plot(res2$nb, type = 'l')
plot(res3$nb, type = 'l')



library(fpopw)
n = 10^6
beta = 2*log(n)

y <- dataGenerator_1D(chpts = n, parameters = 0, type = "gauss")
res1 <- system.time(dust.partitioner.1D(model = "gauss", method = "detIndex_Eval1")$quick(data = y, penalty = beta/2))
res3 <- system.time(fpopw::Fpop(y, beta))
(res1[[1]] - res3[[1]])/res3[[1]]
res1
res3



n = 10^6
beta = 2*log(n)
y <- dataGenerator_1D(chpts = n, parameters = 0, type = "gauss")
res1 <- dust.partitioner.1D(model = "gauss", method = "fastest")$quick(data = y, penalty = beta/2)
res2 <- fpopw::Fpop(y, beta)
res1$nb
res1$changepoints
res2$t.est


########################################################################

n <- 10^8
z <- rnorm(n)
system.time(cs1(z))
system.time(cs2(z))
system.time(cs3(z))
system.time(cs4(z))
system.time(cs5(z))

n <- 10^6
m <- 10
z <- rnorm(n)
system.time(for(i in 1:m){cs1(z)})
system.time(for(i in 1:m){cs2(z)})
system.time(for(i in 1:m){cs3(z)})
system.time(for(i in 1:m){cs4(z)})
system.time(for(i in 1:m){cs5(z)})



n <- 10^5
beta <- 2*log(n)/4
y2 <- dataGenerator_Reg()
y <- y2$x

r1 <- dust.partitioner.meanVar(method = "randIndex_detEval")$quick(data = y, penalty = beta)
r2 <- dust.partitioner.reg(method = "randIndex_randEval")$quick(data = y2, penalty = beta)
r3 <- dust_R_2param_regression(data = y2, penalty = beta, pruningOpt = 2)
r1$changepoints
r2$changepoints
r3$changepoints

r2$costQ - r3$costQ


################################################


n <- 10^4
beta <- 2*log(n)
y <- rnorm(n)
system.time(dust.partitioner.meanVar(method = "fastest")$quick(data = y, penalty = beta))
system.time(dust.partitioner.meanVar(method = "randIndex_randEval")$quick(data = y, penalty = beta))


res1 <- dust.partitioner.1D(model = "gauss", method = "randIndex_detEval")$quick(data = y, penalty = beta)
res2 <- dust.partitioner.meanVar(method = "randIndex_detEval")$quick(data = y, penalty = 2*beta)
res1$lastIndexSet
res2$lastIndexSet
length(res1$lastIndexSet)/n
length(res2$lastIndexSet)/n





################################################################################
################################################################################

n <- 3*10^4
beta <- 2*log(n)
#y <- c(rnorm(n), rnorm(n,sd = 2), rnorm(n, mean = 1, sd = 2))
y <- rnorm(n)
z <- dataGenerator_1D(n, parameters = 40, type = "poisson")
#plot(y, type = 'b')
system.time(dust.partitioner.1D(model = "gauss", method = "fastest")$quick(data = y, penalty = beta))
system.time(dust.partitioner.meanVar(method = "randIndex_randEval")$quick(data = y, penalty = beta))


res1 <- dust.partitioner.1D(model = "gauss", method = "randIndex_detEval")$quick(data = y, penalty = beta)
res2 <- dust.partitioner.meanVar(method = "randIndex_detEval")$quick(data = y, penalty = 2*beta)
res1$lastIndexSet
res2$lastIndexSet
length(res1$lastIndexSet)/n
length(res2$lastIndexSet)/n

res4 <- dust.partitioner.reg(method = "randIndex_detEval")$quick(data = y, penalty = beta)
res3 <- dust_R_2param(y, penalty = beta, pruningOpt = 2)

#all(res2$changepoints == res3$changepoints)
res2$changepoints
res2$lastIndexSet
length(res2$lastIndexSet)
length(res2$lastIndexSet)/(n) * 100

#(res2$costQ - res3$costQ)[n]

################################################################################
################################################################################



res1$lastIndexSet
res2$lastIndexSet


n <- 10^4
beta <- 2 * log(n)/10

type = "poisson"
param = 10
myData <- dataGenerator_1D(n, parameters = param, type = type)
d <- dust.partitioner.1D(model = type, method = "randIndex_randEval", nbLoops = 1)
res <- d$quick(data = myData, penalty = beta)

res2 <- dust_R_1D(myData, type = "poisson", penalty = beta)
res$lastIndexSet
res2$lastIndexSet
length(res$lastIndexSet)/n * 100
length(res2$lastIndexSet)/n * 100


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


test <- function(n = 10^4, type = "gauss", param = 0.5, beta = 0, nbLoops = 5)
{
  myData <- dataGenerator_1D(n, parameters = param, type = type)

  if(type == "binom"){myData <- myData/max(myData)}
  t1 <- system.time(dust.partitioner.1D(model = type, method = "randIndex_randEval", nbLoops = nbLoops)$quick(data = myData, penalty = beta))
  t2 <- system.time(dust.partitioner.1D(model = type, method = "randIndex_detEval", nbLoops = nbLoops)$quick(data = myData, penalty = beta))
  t3 <- system.time(dust.partitioner.1D(model = type, method = "fastest", nbLoops = nbLoops)$quick(data = myData, penalty = beta))
  return(c(t1[[1]],t2[[1]],t3[[1]]))
}


n <- 10^5
beta <- 2 * log(n)
nb <- 10

test(n, "gauss", beta = beta, nbLoops = nb)
test(n, "poisson", beta = beta, nbLoops = nb)
test(n, "exp", beta = beta, nbLoops = nb)
test(n, "geom", beta = beta, nbLoops = nb)
test(n, "bern", beta = beta, nbLoops = nb)
test(n, "binom", beta = beta, nbLoops = nb)
test(n, "negbin", beta = beta, nbLoops = nb)


################################################################################
################################################################################

library(gfpop)
n <- 10^6
beta <- 2 * log(n)
chgtpt <- c(1)
myData <- dataGenerator_1D(n, chgtpt, 0.5, type = "exp")
GRAPH <- graph(type = "std", penalty = beta)

system.time(gfpop(data = myData, type = "exp", mygraph = GRAPH))
system.time(dust.partitioner.1D(model = "exp", method = "fastest")$quick(data = myData, penalty = beta))







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
