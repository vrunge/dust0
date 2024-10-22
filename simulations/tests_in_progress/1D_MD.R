
####

set.seed(45)
d <- 20
data <- dataGenerator_1D(chpts = c(50,100), parameters = c(10,4))
dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)

r1 <- dust.1D(data, penalty = 0.1*log(100))
rM <- dust.MD(dataMD, penalty =  0.1*d*log(100))
rM$costQ/r1$costQ
all(r1$changepoints == rM$changepoints)


d <- 5
data <- dataGenerator_1D()
dataMD <- matrix(rep(data,5), nrow = d, length(data), byrow = T)

r1 <- dust.1D(data)
rM <- dust.MD(dataMD)
r1$costQ/rM$costQ

r1$changepoints
r1$lastIndexSet
r1$nb
r1$costQ

rM$changepoints
rM$lastIndexSet
rM$nb
rM$costQ

####


set.seed(9192)
d <- 5
data <- dataGenerator_1D(chpts = c(50,100), parameters = c(0,1))
dataMD <- matrix(rep(data,5), nrow = d, length(data), byrow = T)

r1 <- dust.1D(data)
rM <- dust.MD(dataMD)

r1$changepoints
rM$changepoints


################################################
########################################################
################################################

d <- 10
data <- dataGenerator_1D(chpts = c(40,60,100), parameters = c(1,0.01,1), type = "poisson")
dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)

r1 <- dust.1D(data, penalty = 0.1*log(100), model = "poisson")
rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), model = "poisson", method = "detIndex_Eval0")
rM$costQ/r1$costQ
r1$changepoints
rM$changepoints

#competitions =>erreur in MD Poisson
segmentation_Cost_1D(data, chpts = r1$changepoints, model = "poisson")
segmentation_Cost_1D(data, chpts = rM$changepoints, model = "poisson")
data

################################################
########################################################
################################################


data <- dataGenerator_1D(chpts = c(50,100), parameters = c(0,1))
data
data_normalization(data, type = "gauss")


d <- 5
data <- dataGenerator_1D(chpts = c(50,100), parameters = c(0,1))
dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)
r1 <- dust.1D(data, penalty = 0.1*log(100))
rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), method = "detIndex_Eval0")
rM$costQ/r1$costQ
r1$changepoints
rM$changepoints


data <- dataGenerator_1D(chpts = c(50,100), parameters = c(10,5), type = "poisson")
dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)
r1 <- dust.1D(data, penalty = 0.1*log(100), model = "poisson")
rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), model = "poisson", method = "detIndex_Eval0")
rM$costQ/r1$costQ
r1$changepoints
rM$changepoints


data <- dataGenerator_1D(chpts = c(50,100), parameters = c(4,1), type = "exp")
dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)
r1 <- dust.1D(data, penalty = 0.1*log(100), model = "exp")
rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), model = "exp", method = "detIndex_Eval0")
rM$costQ/r1$costQ
r1$changepoints
rM$changepoints


data <- dataGenerator_1D(chpts = c(500,1000), parameters = c(0.7,0.6), type = "geom")
dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)
r1 <- dust.1D(data, penalty = 1*log(1000), model = "geom")
rM <- dust.MD(dataMD, penalty =  1*d*log(1000), model = "geom", method = "detIndex_Eval0")
rM$costQ/r1$costQ
r1$changepoints
rM$changepoints



data <- dataGenerator_1D(chpts = c(50,100), parameters = c(0.7,0.6), type = "binom")
data <- data_normalization(data, type = "binom")
dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)
r1 <- dust.1D(data, penalty = 0.1*log(100), model = "binom")
rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), model = "binom", method = "detIndex_Eval0")
rM$costQ/r1$costQ
r1$changepoints
rM$changepoints













### OK

data <- dataGenerator_1D(chpts = c(50,100), parameters = c(0.7,0.6), type = "bern")
dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)
r1 <- dust.1D(data, penalty = 0.1*log(100), model = "bern")
rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), model = "bern")
rM$costQ/r1$costQ
r1$changepoints
rM$changepoints


### PAS OK A CAUSE DE method = "detIndex_Eval0"

data <- dataGenerator_1D(chpts = c(50,100), parameters = c(0.7,0.6), type = "bern")
dataMD <- matrix(rep(data, d), nrow = d, length(data), byrow = T)
r1 <- dust.1D(data, penalty = 0.1*log(100), model = "bern")
rM <- dust.MD(dataMD, penalty =  0.1*d*log(100), model = "bern", method = "detIndex_Eval1")
rM$costQ/r1$costQ
r1$changepoints
rM$changepoints







