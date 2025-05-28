

getwd()





library(dust)
###
### GAUSS 3D
###

n <- 1000
data <- dataGenerator_MD(chpts = n,
                         parameters = data.frame(ts1 = 0, ts2 = 0, ts3 = 0))
res <- dust.MD(data, method = "detIndex_Eval6", model = "gauss", nbLoops = 100)






