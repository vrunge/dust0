
n = 1e2
beta = 2*log(n)
y <- rnorm(n)
dust.partitioner.1D()$quick(data = y, penalty = beta)




res <- dust.partitioner.1D()$quick(rnorm(1e4))
res$changepoints
res$lastIndexSet


DUST = new(DUST_1D, "gauss", "fastest", NULL) # max 1D dual Gauss, indice last one (time 1)
DUST = new(DUST_1D, "gauss", "deterministic", NULL) # idem
DUST = new(DUST_1D, "gauss", "full.random", NULL) # random in indices and dual random eval unif (0,1) (time 3)
DUST = new(DUST_1D, "gauss", "half.random", NULL) # random in indices and eval dual max (time 2)


DUST$init_raw(rnorm(1e3), NULL)
DUST$compute()
DUST$get_partition()

res <- dust.partitioner.1D("gauss", "full.random")$quick(c(rnorm(1e3), rnorm(1e3, mean = 1)))


#################################################################


n = 1e2
beta = 2 * log(n)
y <- rnorm(n)

system.time(dust.partitioner.1D()$quick(data = y, penalty = beta))
system.time(dust.partitioner.1D("gauss", "half.random")$quick(data = y, penalty = beta))
system.time(fpopw::Fpop(y, beta))



#################################################################

library(microbenchmark)
library(fpopw)

n = 1e4
K = 100
Y = matrix(rnorm(n * K), ncol = K)
dim(Y)
beta = 2 * log(n)

j = 1
jj = 0
microbenchmark::microbenchmark(
  DUST = dust.partitioner.1D()$quick(y, beta),
  FPOP = fpopw::Fpop(y, beta),
  times = K,
  setup = {
    y = Y[, j]
    j = j + jj %% 2
    jj = j + 1
  }
)



#################################################################
### POISSON ###

n = 10^3
beta = 2*log(n)/2
y <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(1,2), type = "poisson")
plot(y)
res1 <- dust.partitioner.1D(model = "poisson", method = "all.random")$quick(data = y, penalty = beta)
res2 <- dust_R_1D(y, type = "poisson", penalty = beta, pruningOpt = 0)


y

res1$changepoints
res2$changepoints
length(res1$lastIndexSet)

(res1$costQ - res2$costQ)[n]

#################################################################
### EXP ###

n = 1e3
beta = 2*log(n)/10
y <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(0.1,200), type = "exp")
plot(y)
res1 <- dust.partitioner.1D(model = "exp", method = "all.random")$quick(data = y, penalty = beta)
res2 <- dust_R_1D(y, type = "exp", penalty = beta, pruningOpt = 0)

res1$changepoints
res2$changepoints

length(res1$lastIndexSet)
(res1$costQ - res2$costQ)[n]



#################################################################
### GEOM ###

n = 1e3
beta = log(n)/1
y <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(0.86,0.92), type = "geom")
y;plot(y)
res1 <- dust.partitioner.1D(model = "geom", method = "all.random")$quick(data = y, penalty = beta)
res2 <- dust_R_1D(y, type = "geom", penalty = beta, pruningOpt = 0)

res1$changepoints
res2$changepoints

length(res1$lastIndexSet)
(res1$costQ - res2$costQ)[n]




#################################################################
### BERN ###

n = 10^3
beta = log(n)
y <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(0.1,0.5), type = "bern")
y; plot(y)
res1 <- dust.partitioner.1D(model = "bern", method = "all.random")$quick(data = y, penalty = beta)
res2 <- dust_R_1D(y, type = "bern", penalty = beta, pruningOpt = 0)



res1$lastIndexSet
res1$changepoints
res2$changepoints

(res1$costQ - res2$costQ)[n]
length(res1$lastIndexSet)





