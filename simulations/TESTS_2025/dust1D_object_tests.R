

#################################################################################
#################################################################################
### TEST dust.1D


model = "gauss"
method = "randIndex_Eval2"
nbLoops = 10


res <- dust.1D(data = data2, penalty = penalty)
res$changepoints
plot(data2)
plot(res$costQ)

######


obj_dust <- new(DUST_1D, "variance", "randIndex_Eval4", 5)
penalty <- 2*log(5)
data_all <- NULL
for(i in 1:5)
{
  data <- dataGenerator_1D(chpts = c(500,1000), parameters = c(1,1.5), type = "variance")
  data_all <- c(data_all, data)
  obj_dust$append_c(data, penalty)
  obj_dust$update_partition()
}
resObject <- obj_dust$get_partition()
res <- dust.1D(data = data_all,
               penalty = penalty,
               model = "variance",
               method = "randIndex_Eval4", nbLoops = 5)
all(res$changepoints == resObject$changepoints)
all(res$costQ == resObject$costQ)

res$changepoints
plot(data_all)



obj_dust$append_c(data, penalty)
obj_dust$update_partition()
res2 <- obj_dust$get_partition()
res2$changepoints
plot(res2$costQ)

all(res2$costQ == res$costQ)

res$changepoints


####################################





