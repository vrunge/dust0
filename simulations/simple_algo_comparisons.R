
###############################
############ GAUSS ############
###############################
#devtools::install_github("vrunge/gfpop")
#devtools::install_github("vrunge/dust")
library(gfpop)
library(fpopw)
library(dust)

###############################
###############################

n <- 10^6
data <- dataGenerator_1D(n, parameters = 0, sdNoise = 1, type = "gauss")
#plot(data)
pen <- 2*log(n)/2
res_fpop <- fpopw::Fpop(data, lambda = pen)
res_fpop$t.est

res_dust <- dust.1D(data, penalty = pen/2, model = "gauss", method = "detIndex_Eval1")
res_dust$changepoints

all(res_dust$changepoints == res_fpop$t.est)


system.time(res_fpop <- fpopw::Fpop(data, lambda = pen))
system.time(res_dust <- dust.1D(data, penalty = pen/2, model = "gauss", method = "detIndex_Eval1"))


#################################
############ POISSON ############
#################################

n <- 10^5
chgtpt <- c(0.1, 0.3, 0.5, 0.8, 1)
data <- dataGenerator_1D(chpts = chgtpt*n, parameters = c(5, 1, 5, 1, 5), type = "poisson")
plot(data)

pen <- 2 * log(n)
GRAPH <- graph(type = "std", penalty = pen)
res_GFPOP <- gfpop(data = data, type = "poisson", mygraph = GRAPH)
res_DUST <- dust.1D(data, penalty = pen, model = "poisson", method = "detIndex_Eval4")

res_GFPOP$globalCost + pen*(length(chgtpt)-1)
res_DUST$costQ[n]

res_GFPOP$changepoints
res_DUST$changepoints

system.time(gfpop(data = data, type = "poisson", mygraph = GRAPH))
system.time(dust.1D(data, penalty = pen, model = "poisson", method = "detIndex_Eval4"))



################################
############ NEGBIN ############
################################

library()

n <- 10^6
chgtpt <- c(0.1, 0.3, 0.5, 0.8, 1)
data <- dataGenerator_1D(chpts = chgtpt*n, parameters = c(5, 1, 5, 1, 5), type = "poisson")
#plot(data)

pen <- 2 * log(n)
GRAPH <- graph(type = "std", penalty = pen)
res_GFPOP <- gfpop(data = data, type = "poisson", mygraph = GRAPH)
res_DUST <- dust.1D(data, penalty = pen, model = "poisson", method = "detIndex_Eval4")

res_GFPOP$globalCost + pen*(length(chgtpt)-1)
res_DUST$costQ[n]

res_GFPOP$changepoints
res_DUST$changepoints


####################################
####################################

n <- 10^7
data <- dataGenerator_1D(chpts = n, parameters = c(2), type = "poisson")
system.time(gfpop(data = data, type = "poisson", mygraph = GRAPH))
system.time(dust.1D(data, penalty = pen, model = "poisson", method = "detIndex_Eval4"))



########################################################################
####################################

n <- 10^4
chgtpt <- c(0.1, 0.3, 0.5, 0.8, 1)
data <- dataGenerator_1D(chpts = chgtpt*n, parameters = c(5, 1, 5, 1, 5), type = "poisson")
plot(data)

pen <- 2 * log(n)
GRAPH <- graph(type = "std", penalty = pen)
res_GFPOP <- gfpop(data = data, type = "poisson", mygraph = GRAPH)
res_DUST4 <- dust.1D(data, penalty = pen, model = "poisson", method = "detIndex_Eval4")
res_DUST1 <- dust.1D(data, penalty = pen, model = "poisson", method = "detIndex_Eval1")

res_GFPOP$globalCost + pen*(length(chgtpt)-1)
res_DUST4$costQ[n]
res_DUST1$costQ[n]
plot(res_DUST4$nb)
plot(res_DUST1$nb)

plot(res_DUST4$nb - res_DUST1$nb)

res_GFPOP$changepoints
res_DUST4$changepoints
res_DUST1$changepoints


n <- 10^6
chgtpt <- 1
data <- dataGenerator_1D(chpts = chgtpt*n, parameters = 5, type = "poisson")
#plot(data)
system.time(res_GFPOP <- gfpop(data = data, type = "poisson", mygraph = GRAPH))
system.time(res_DUST4 <- dust.1D(data, penalty = pen, model = "poisson", method = "detIndex_Eval4"))
system.time(res_DUST1 <- dust.1D(data, penalty = pen, model = "poisson", method = "detIndex_Eval1"))

res_DUST4$nb[length(res_DUST4$nb)]
res_DUST1$nb[length(res_DUST1$nb)]

plot(res_DUST4$nb)




