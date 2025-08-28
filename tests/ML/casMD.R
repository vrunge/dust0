

getwd()





library(dust0)
###
### GAUSS 3D
###

n <- 4000
data <- dataGenerator_MD(chpts = n,
                         parameters = data.frame(ts1 = 0, ts2 = 0, ts3 = 0))
res <- dust.MD(data, method = "randIndex_Eval6", model = "gauss", constraints_l = 2, nbLoops = 100)
plot(res$nb)





