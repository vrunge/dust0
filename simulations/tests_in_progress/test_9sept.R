library(dust0)

tpe <- "exp"
n <- 10^4
beta <- 2*log(n)
y <- dataGenerator_1D(c(n/2, n), parameters = c(0.2, 0.3), type = tpe)
res1 <- dust.partitioner.1D(method = "detIndex_Eval4", model = tpe)$quick(data = y, penalty = beta)
res2 <- dust.partitioner.1D(method = "detIndex_Eval2", model = tpe)$quick(data = y, penalty = beta)


a <- dust.partitioner.1D(method = "detIndex_Eval4", model = tpe)
a$init(y, beta)
a$compute(y)
res <- a$get_partition()
res$changepoints
####

?DUST_1D
?dust.partitioner.1D


########################################

n <- 10^2
beta <- 2*log(n)
y <- dataGenerator_1D(n, parameters = 0, type = "gauss")
a <- dust.partitioner.1D(method = "randIndex_Eval5", model = "gauss")$quick(data = y, penalty = beta)
a$changepoints
a$costQ
length(a$costQ)

a$lastIndexSet
length(a$lastIndexSet)

a$nb
length(a$nb)


