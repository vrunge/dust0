
tpe <- "poisson"
n <- 10^5
beta <- 2*log(n)/5
y <- dataGenerator_1D(n, parameters = 0, type = tpe)
res1 <- dust.partitioner.1D(method = "detIndex_Eval4", model = tpe)$quick(data = y, penalty = beta)
res2 <- dust.partitioner.1D(method = "detIndex_Eval2", model = tpe)$quick(data = y, penalty = beta)

all(res1$changepoints == res2$changepoints)
res2$changepoints
plot(res1$nb - res2$nb)

n <- 10^4
beta <- 2*log(n)
y <- dataGenerator_1D(n, parameters = 10, type = tpe)
system.time(dust.partitioner.1D(method = "randIndex_Eval2", model = tpe)$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval4", model = tpe)$quick(data = y, penalty = beta))


n <- 10^3
beta <- 2*log(n)
y <- dataGenerator_1D(n, parameters = 0, type = "gauss")
a <- dust.partitioner.1D(method = "randIndex_Eval2", model = "gauss")$quick(data = y, penalty = beta)
b <- dust.partitioner.1D(method = "randIndex_Eval4", model = "gauss")$quick(data = y, penalty = beta)

a$changepoints
b$changepoints

a$lastIndexSet
b$lastIndexSet





####################################################################################
####################################################################################
####################################################################################

###
###
### testing all methods on Gaussian DATA
###
###
n <- 10^5
beta <- 2*log(n)
y <- rnorm(n)
res0 <- dust.partitioner.1D(method = "randIndex_Eval0")$quick(data = y, penalty = beta)
res1 <- dust.partitioner.1D(method = "randIndex_Eval1")$quick(data = y, penalty = beta)
res2 <- dust.partitioner.1D(method = "randIndex_Eval2")$quick(data = y, penalty = beta)
res3 <- dust.partitioner.1D(method = "randIndex_Eval3")$quick(data = y, penalty = beta)
res4 <- dust.partitioner.1D(method = "randIndex_Eval4")$quick(data = y, penalty = beta)
#res5 <- dust.partitioner.1D(method = "randIndex_Eval5")$quick(data = y, penalty = beta)
plot(res2$nb-res5$nb)

res5 <- dust.partitioner.1D(method = "detIndex_Eval0")$quick(data = y, penalty = beta)
res6 <- dust.partitioner.1D(method = "detIndex_Eval1")$quick(data = y, penalty = beta)
res7 <- dust.partitioner.1D(method = "detIndex_Eval2", nbLoops = 4)$quick(data = y, penalty = beta)

length(res1$lastIndexSet)
length(res2$lastIndexSet)
length(res3$lastIndexSet)

length(res5$lastIndexSet)
length(res6$lastIndexSet)
length(res7$lastIndexSet)

res1$nb[n]/n * 100
res2$nb[n]/n * 100

(n*(n+1)/2) /sum(res1$nb)
(n*(n+1)/2) /sum(res2$nb)

sum(res1$nb)/sum(res2$nb)



###
###
### testing all methods on Poisson DATA
###
###
n <- 10^5
beta <- 4*log(n)/6
y <- dataGenerator_1D(n, parameters = 4, type = "poisson")
res0 <- dust.partitioner.1D(method = "randIndex_Eval0", model = "poisson")$quick(data = y, penalty = beta)
res2 <- dust.partitioner.1D(method = "randIndex_Eval2", model = "poisson")$quick(data = y, penalty = beta)
res4 <- dust.partitioner.1D(method = "randIndex_Eval4", model = "poisson")$quick(data = y, penalty = beta)

res5 <- dust.partitioner.1D(method = "detIndex_Eval0", model = "poisson")$quick(data = y, penalty = beta)
res7 <- dust.partitioner.1D(method = "detIndex_Eval2", model = "poisson")$quick(data = y, penalty = beta)
res9 <- dust.partitioner.1D(method = "detIndex_Eval4", model = "poisson")$quick(data = y, penalty = beta)

length(res0$lastIndexSet)
length(res2$lastIndexSet)
length(res4$lastIndexSet)

length(res5$lastIndexSet)
length(res7$lastIndexSet)
length(res9$lastIndexSet)

sum(res7$costQ - res0$costQ)


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
#system.time(dust.partitioner.1D(method = "randIndex_Eval0")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval1")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval2")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval3")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval4")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval5")$quick(data = y, penalty = beta))


#system.time(dust.partitioner.1D(method = "detIndex_Eval0")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval1")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval2")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval3")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval4")$quick(data = y, penalty = beta))

###
###
### TIME all methods on Poisson DATA
###
###
n <- 10^6
beta <- 2*log(n)
y <- dataGenerator_1D(n, parameters = 50, type = "poisson")
#system.time(dust.partitioner.1D(method = "randIndex_Eval0", model = "poisson")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval2", model = "poisson")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval3", model = "poisson")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "randIndex_Eval4", model = "poisson")$quick(data = y, penalty = beta))

#system.time(dust.partitioner.1D(method = "detIndex_Eval0", model = "poisson")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval2", model = "poisson")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval3", model = "poisson")$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D(method = "detIndex_Eval4", model = "poisson")$quick(data = y, penalty = beta))


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
r1$nb[n]/n * 100
r2$nb[n]/n * 100

####################################################################################
####################################################################################
###
###
### meanvar
###
###

n <- 10^4
beta <- 2*log(n)
y <- c(dataGenerator_meanVar(chpts = c(n/4,n/2,n), means = c(0,1,0), sds = c(1,1,1.2)),rep(0, n))
#y <- rep(0, n)
plot(y)
r1 <- dust.partitioner.meanVar(method = "randIndex_Eval0")$quick(data = y, penalty = beta)
r2 <- dust.partitioner.meanVar(method = "randIndex_Eval2")$quick(data = y, penalty = beta)
r3 <- dust.partitioner.meanVar(method = "randIndex_Eval4")$quick(data = y, penalty = beta)
r4 <- dust.partitioner.meanVar(method = "detIndex_Eval0")$quick(data = y, penalty = beta)
r5 <- dust.partitioner.meanVar(method = "detIndex_Eval2")$quick(data = y, penalty = beta)
r6 <- dust.partitioner.meanVar(method = "detIndex_Eval4")$quick(data = y, penalty = beta)
r1$changepoints
r2$changepoints
r3$changepoints
r4$changepoints
r5$changepoints
r6$changepoints
r1$nb[n]
r2$nb[n]
r3$nb[n]
r4$nb[n]
r5$nb[n]
r6$nb[n]
r1$nb[n]/n * 100
r2$nb[n]/n * 100
r3$nb[n]/n * 100
r4$nb[n]/n * 100
r5$nb[n]/n * 100
r6$nb[n]/n * 100




n <- 10^3
beta <- 2*log(n)
y <-   dataGenerator_meanVar(chpts = n)
plot(y)
system.time(dust.partitioner.meanVar(method = "randIndex_Eval0")$quick(data = y, penalty = beta))
system.time(dust.partitioner.meanVar(method = "randIndex_Eval2")$quick(data = y, penalty = beta))
system.time(dust.partitioner.meanVar(method = "randIndex_Eval4")$quick(data = y, penalty = beta))



