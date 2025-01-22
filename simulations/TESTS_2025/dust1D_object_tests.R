

#################################################################################
#################################################################################
### TEST dust.1D


data = dataGenerator_1D(chpts = c(500,1000), parameters = c(0,1), sdNoise = 0.2, type = "gauss")

penalty = 2*log(length(data))/200
data2 <- rep(data, 2)
model = "gauss"
method = "randIndex_Eval2"
nbLoops = 10


res <- dust.1D(data = data2, penalty = penalty)
res$changepoints
plot(data2)
plot(res$costQ)

######
obj_dust <- new(DUST_1D, model, method, nbLoops)

obj_dust$append_c(data, penalty)
obj_dust$update_partition()

obj_dust$append_c(data, penalty)
obj_dust$update_partition()
res2 <- obj_dust$get_partition()
res2$changepoints
plot(res2$costQ)

all(res2$costQ == res$costQ)

res$changepoints


####################################





