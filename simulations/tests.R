

#dataGenerator

dataGenerator(chpts = c(50,100), parameter = c(0,1), type = "gauss")
plot(dataGenerator(chpts = c(50,100), parameter = c(0,1), type = "gauss"))


plot(dataGenerator(chpts = c(30,100,120), parameter = c(1,10,1), type = "poisson"))
data <- dataGenerator(chpts = c(300,1000,1200), parameter = c(1,1.4,1), type = "gauss")

type <- "poisson"
s1 <- 4
s2 <- 1
t <- 10
x <- dataGenerator_1D(chpts = c(15,30), parameters = c(1,10), type = type)

mu = seq(0,1,length.out = 100)
res <- dual_1D(mu = mu, x = x, s1 = s1, s2 = s2, t = t, type = type)
plot(mu, res, type = 'l')

res <- dual_1D(mu = mu, x = x, s1 = s1, s2 = s2, t = t, type = type, penalty = 2*log(30), OP = TRUE)
plot(mu, res, type = 'l')

OP_R(x, penalty = 2*log(30), type = type)


x <- dataGenerator_1D(chpts = c(15,30), parameters = c(1,10), type = "gauss")
dual_1D(mu =  1:99/100, x = x, s1 = 4, s2 = 2, t = 10, type = "gauss")


#####################
##### dust_R_1D #####
#####################

n <- 100
type <- "gauss"
data <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(0,2), type = type)
dust_R_1D(data, 2*log(n), type = type)

type <- "poisson"
data <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(5,2), type = type)
dust_R_1D(data, 2*log(n), type = type, pruningOpt = 1)


type <- "exp"
data <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(5,1), type = type)
dust_R_1D(data, 2*log(n), type = type)

type <- "bern"
n <- 100
data <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(0.4,0.6), type = type)
plot(data)
dust_R_1D(data, log(n)/2, type = type)




diff(res$lastIndexSet[-1])
res$changepoints
res2 <- OP_R_1D(data, 2*log(n), type = type)
res2$changepoints

res$costQ == res2$costQ

#####################
n <- 1000
type <- "poisson"
data <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(5,1), type = type)
set.seed(1)
res <- dust_R_1D(data, 2*log(n), type = type, pruningOpt = 1)
res$changepoints
plot(data, type = 'l')
diff(res$lastIndexSet[-1])
res$lastIndexSet
res2 <- OP_R_1D(data, 2*log(n), type = type)
res2$changepoints

res$costQ == res2$costQ

#########################################################################################################


n <- 1000
data <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(0.4,0.5), type = "geom")
res <- dust_R_1D(data, 5*log(n), type = "geom")
res$changepoints
plot(data, type = 'l')
diff(res$lastIndexSet[-1])
res$lastIndexSet




n <- 10000
data <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(4,6), type = "exp")
res <- dust_R_1D(data, 1*log(n))
res$changepoints
plot(data, type = 'l')
diff(res$lastIndexSet[-1])
res$lastIndexSet


n <- 1000
data <- dataGenerator_1D(chpts = c(n/2,n), parameters = c(0.3,0.5), type = "bern")
res <- dust_R_1D(data, log(n))
res$changepoints
plot(data, type = 'l')
diff(res$lastIndexSet[-1])
res$lastIndexSet




#### VERY BIG TEST

n <- 10^5
type <- "poisson"
data <- dataGenerator_1D(chpts = n, parameter = 5, type = type)
a <- Sys.time()
res <- dust_R_1D(data, 4*log(n), type = type, pruningOpt = 2)
b <- Sys.time()
b-a
res$changepoint
res$lastIndexSet
length(res$lastIndexSet)/n*100



nb <- 10^4
s1 <- 18
s2 <- 15
j <- 1
bar <- rep(0,99)
for(i in 1:nb)
{
    data <- dataGenerator_1D(chpts = 20, parameters = 0, type = "gauss")
    res <- dual_1D(mu =  1:99/100, data = data, s1 = s1, s2 = s2, t = 20,
                   type = "gauss", OP = TRUE, penalty = 1)
    vertical <- which(res$dualValues > res$pruningBound)


    if(!(s1 %in% res$lastIndexSet)){j<-j+1}
    if(s1 %in% res$lastIndexSet)
    {
    u <- rep(1,99)
    bar[vertical] <- bar[vertical] + u[vertical]
    }
}
print(j)
barplot(bar)


##############

data <- dataGenerator_1D(chpts = 10, parameters = 0.5, type = "bern")
data
plot_dual_1D(mu = 1:99/100, data = data, s1 = 4, s2 = 1, type = "bern")



data <- dataGenerator_1D(chpts = 10, parameters = 0.8, type = "geom")
data
plot_dual_1D(mu = 1:999/1000, data = data, s1 = 4, s2 = 3, type = "geom")




data <- dataGenerator_1D(chpts = 10, parameters = 4, type = "poisson")
data
plot_dual_1D(mu = 1:999/1000, data = data, s1 = 4, s2 = 1, type = "poisson")


data <- dataGenerator_1D(chpts = 10, parameters = 1, type = "exp")
data
plot_dual_1D(mu = 1:999/1000, data = data, s1 = 4, s2 = 3, type = "exp")


data <- dataGenerator_1D(chpts = 10, parameters = 1, type = "gauss")
data
plot_dual_1D(mu = 1:99/100, data = data, s1 = 4, s2 = 3, type = "gauss")



#####################################################################################
#####################################################################################
#####################################################################################


p <- 2
data <- dataGenerator_1D(100, sdNoise = 2) # noise != 1 to generate changes at random positions
pen <- 2*log(100)
op1D <- OP_R_1D(data = data, penalty = pen, type = "gauss")
dataM <- matrix(data, p, 100, byrow = T)
opMD <- OP_R_MultiD(data = dataM, penalty = p*pen, type = "gauss")

opMD2 <- dust_R_MultiD(data = dataM, penalty = p*pen, type = "gauss", pruningOpt = 2)


op1D$changepoints
opMD$changepoints
opMD2$changepoints
opMD2$nb

abs(op1D$costQ - opMD$costQ/p) < 10^(-13)
abs(op1D$costQ - opMD2$costQ/p) < 10^(-13)




p <- 20
data <- dataGenerator_1D(100, parameter = 50, type = "poisson") # noise != 1 to generate changes at random positions
pen <- 1
op1D <- OP_R_1D(data = data, penalty = pen, type = "poisson")

dataM <- matrix(data, p, 100, byrow = T)
opMD <- OP_R_MultiD(data = dataM, penalty = p*pen, type = "poisson")
opMD2 <- dust_R_MultiD(data = dataM, penalty = p*pen, type = "poisson", pruningOpt = 2)


op1D$changepoints
opMD$changepoints
opMD2$changepoints
opMD2$nb

abs((op1D$costQ - opMD$costQ/p)/op1D$costQ) < 10^(-13)
abs((op1D$costQ - opMD2$costQ/p)/op1D$costQ) < 10^(-13)

