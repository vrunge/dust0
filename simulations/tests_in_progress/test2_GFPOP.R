

#devtools::install_github("vrunge/gfpop")
library(gfpop)


n <- 10^5
beta <- 2 * log(n)
chgtpt <- c(1)
myData <- dataGenerator(n, 1, 50, type = "poisson")
myData <- dataGenerator_1D(n, parameters = 0.5, type = "binom")
#plot(myData)

GRAPH <- graph(type = "std", penalty = beta)

#system.time(gfpop(data = myData, type = "binom", mygraph = GRAPH))
system.time(dust.partitioner.1D(model = "binom", method = "randIndex_randEval")$quick(data = myData, penalty = beta))
system.time(dust.partitioner.1D(model = "binom", method = "randIndex_detEval")$quick(data = myData, penalty = beta))


##########################################
##########################################


n <- 1000
beta <- 2 * log(n)
chgtpt <- c(0.1, 0.3, 0.5, 0.8, 1)
myData <- dataGenerator(n, chgtpt, c(10, 1, 10, 1, 10), type = "poisson")
#plot(myData)

GRAPH <- graph(type = "std", penalty = beta)

resGFPOP <- gfpop(data = myData, type = "poisson", mygraph = GRAPH)
resDUST <- dust.partitioner.1D(model = "poisson", method = "randIndex_randEval")$quick(data = myData, penalty = beta)


resGFPOP$globalCost
resDUST$costQ[n]


resGFPOP$changepoints
resDUST$changepoints

##########################################
##########################################
##########################################
##########################################
##########################################

n <- 10^5
beta <- 2 * log(n)
chgtpt <- c(1)
myData <- dataGenerator(n, 1, 50, type = "poisson")
#plot(myData)

GRAPH <- graph(type = "std", penalty = beta)

system.time(gfpop(data = myData, type = "poisson", mygraph = GRAPH))
system.time(dust.partitioner.1D(model = "poisson", method = "randIndex_randEval")$quick(data = myData, penalty = beta))
system.time(dust.partitioner.1D(model = "poisson", method = "randIndex_detEval")$quick(data = myData, penalty = beta))


### ### ### ### ### ### ###

gfpop(data = myData, type = "poisson", mygraph = GRAPH)



