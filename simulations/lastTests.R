

###
###
### testing all methods on Gaussian DATA
###
###
n <- 10^5
beta <- 2*log(n)
y <- rnorm(n)
res1 <- dust.partitioner.1D(method = "randIndex_Eval0")$quick(data = y, penalty = beta)
res2 <- dust.partitioner.1D(method = "randIndex_Eval1")$quick(data = y, penalty = beta)
res3 <- dust.partitioner.1D(method = "randIndex_Eval2", nbLoops = 4)$quick(data = y, penalty = beta)

res5 <- dust.partitioner.1D(method = "detIndex_Eval0")$quick(data = y, penalty = beta)
res6 <- dust.partitioner.1D(method = "detIndex_Eval1")$quick(data = y, penalty = beta)
res7 <- dust.partitioner.1D(method = "detIndex_Eval2", nbLoops = 4)$quick(data = y, penalty = beta)

length(res1$lastIndexSet)
length(res2$lastIndexSet)
length(res3$lastIndexSet)

length(res5$lastIndexSet)
length(res6$lastIndexSet)
length(res7$lastIndexSet)



###
###
### testing all methods on Poisson DATA
###
###
n <- 10^5
beta <- 4*log(n)
y <- dataGenerator_1D(n, parameters = 4, type = "poisson")
res1 <- dust.partitioner.1D(method = "randIndex_Eval0", model = "poisson")$quick(data = y, penalty = beta)
res3 <- dust.partitioner.1D(method = "randIndex_Eval2", model = "poisson", nbLoops = 20)$quick(data = y, penalty = beta)

res5 <- dust.partitioner.1D(method = "detIndex_Eval0", model = "poisson")$quick(data = y, penalty = beta)
res7 <- dust.partitioner.1D(method = "detIndex_Eval2", model = "poisson", nbLoops = 20)$quick(data = y, penalty = beta)

length(res1$lastIndexSet)
length(res3$lastIndexSet)

length(res5$lastIndexSet)
length(res7$lastIndexSet)


####################################################################################
####################################################################################

###
###
### TIME all methods on Gaussian DATA
###
###
n <- 10^6
beta <- 2*log(n)
y <- rnorm(n)
system.time(dust.partitioner.1D(method = "randIndex_Eval0")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval1")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval2")$quick(data = y, penalty = beta))

system.time(dust.partitioner.1D(method = "detIndex_Eval0")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval1")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval2")$quick(data = y, penalty = beta))


###
###
### TIME all methods on Poisson DATA
###
###
n <- 10^6
beta <- 2*log(n)
y <- dataGenerator_1D(n, parameters = 4, type = "poisson")
system.time(dust.partitioner.1D(method = "randIndex_Eval0", model = "poisson")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval2", model = "poisson")$quick(data = y, penalty = beta))

system.time(dust.partitioner.1D(method = "detIndex_Eval0", model = "poisson")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval2", model = "poisson")$quick(data = y, penalty = beta))


####################################################################################
####################################################################################
###
###
### regression
###
###

n <- 10^4
beta <- 2*log(n)
y <- dataGenerator_Reg(chpts = c(n/4,n), A = c(2,-1),  B = c(-1,2), meansX = c(10,0.1), type = "simple")
r1 <- dust.partitioner.reg(method = "randIndex_Eval0")$quick(data = y, penalty = beta)
r2 <- dust.partitioner.reg(method = "randIndex_Eval2")$quick(data = y, penalty = beta)
r1$changepoints
r2$changepoints
r1$nb[n]/n
r2$nb[n]/n

####################################################################################
####################################################################################
###
###
### meanvar
###
###

n <- 10^4
beta <- 2*log(n)
y <-   dataGenerator_MV(chpts = c(n/4,n/2,n), means = c(0,1,0), sds = c(1,1,2))
r1 <- dust.partitioner.meanVar(method = "randIndex_Eval0")$quick(data = y, penalty = beta)
r2 <- dust.partitioner.meanVar(method = "randIndex_Eval2")$quick(data = y, penalty = beta)
r1$changepoints
r2$changepoints
r1$nb[n]/n
r2$nb[n]/n





