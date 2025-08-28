
#install.packages("fpop", repos="http://R-Forge.R-project.org")
library(fpop)
for(i in seq(from = 7, to = 9, by = 0.1))
{
  print(i)
  n <- i*10^7
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = c(0,1,0,1,0,1,0,1,0,1),
                           type = "gauss",
                           sdNoise = 1)
  print(system.time(multiBinSeg(data, Kmax = 9)))
}

# response 75.10^6

################################################################################

library(fpopw)
for(i in seq(from = 3, to = 5, by = 0.1))
{
  print(i)
  n <- i*10^7
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = c(0,1,0,1,0,1,0,1,0,1),
                           type = "gauss",
                           sdNoise = 1)
  print(system.time(fpopw::Fpop(x = data, lambda = 2*log(n))))
}

# response 30.10^6

################################################################################

library(dust0)
for(i in seq(from = 4, to = 5, by = 0.1))
{
  print(i)
  n <- i*10^7
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = c(0,1,0,1,0,1,0,1,0,1),
                           type = "gauss",
                           sdNoise = 1)
  print(system.time(dust.gauss.1D(data, penalty = 2*log(n)/2, model = "gauss")))
}

# response 42.10^6

################################################################################

for(i in seq(from = 1, to = 1.5, by = 0.1))
{
  print(i)
  n <- i*10^5
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = c(0,1,0,1,0,1,0,1,0,1),
                           type = "gauss",
                           sdNoise = 1)
  print(system.time(flat_OP_1D(data,inPenalty = 2*log(n)/2)))
}

# response 120 10^3

42*10^6 /(120 *10^3)


################################################################################


n <- 2*10^3
cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
data <- dataGenerator_1D(chpts = cpts,
                         parameters = c(0,1,0,1,0,1,0,1,0,1),
                         type = "gauss",
                         sdNoise = 1)

res1 <- flat_OP_1D(data,inPenalty = 2*log(n)/2)
res1$changepoints
res2 <- dust.gauss.1D(data, penalty = 2*log(n)/2, model = "gauss")
res2$changepoints
res1$changepoints == res2$changepoints

#####################################################################################
################################################################################
################################################################################





################################################################################
library(dust0)
library(fpopw)

n <- 4
cpts <- (1:10)*10^(n-1)
data <- dataGenerator_1D(chpts = cpts,
                         parameters = c(0,1,0,1,0,1,0,1,0,1),
                         type = "gauss",
                         sdNoise = 1)

system.time(dust.gauss.1D(data, penalty = 2*log(10^n)/2, model = "gauss"))
system.time(Fpop(x = data, lambda = 2*log(10^n)))
#system.time(flat_OP_1D(data,inPenalty = 2*log(10^n)/2))

res1 <- Fpop(x = data, lambda = 2*log(10^n))
res1$t.est

################################################################################
library(changepoint)

result <- cpt.mean(data, method = "BinSeg", Q = 9)
cpts(result)
res1$t.est == cpts(result)

################################################################################
library(wbs)
wbs_result <- wbs(data, N = 100)
cpt_model <- changepoints(wbs_result, Kmax = 9)
cpt_model$cpt.th
res1$t.est == c(sort(cpt_model$cpt.th[[1]]), 10^n)

################################################################################
library(not)
# Apply NOT (SBS is a special case)
not_result <- not(data, contrast = "pcwsConstMean", M = 10000)

res1$t.est == c(sort(not_result$solution.path$cpt[[10]]),10^n)




