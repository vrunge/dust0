




res <- dust.partitioner.1D()$quick(c(rnorm(1e3), rnorm(1e3, mean = 1)))
res$changepoints




DUST = new(DUST_1D, "gauss", "fastest", NULL) # max 1D dual Gauss, indice last one (time 1)

DUST = new(DUST_1D, "gauss", "deterministic", NULL) # idem
DUST = new(DUST_1D, "gauss", "full.random", NULL) # random in indices and dual random eval unif (0,1) (time 3)
DUST = new(DUST_1D, "gauss", "half.random", NULL) # random in indices and eval dual max (time 2)


DUST$init_raw(rnorm(1e3), NULL)

DUST$compute()

DUST$get_partition()



res <- dust.partitioner.1D("gauss", "full.random")$quick(c(rnorm(1e3), rnorm(1e3, mean = 1)))

#####

library(microbenchmark)
library(fpopw)


n = 1e6
K = 10
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



