

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

0  497  499  500  976  979  980  985  986  987  988  989  990  992  993  994  995  996  997  998  999 1000
497  499  500  986  987  988  990  993  994  996  997  999 1000
497  499  500  986  987  988  990  993  994  996  997  999 1000
497  499  500  978  985  986  988  989  990  993  994  995  997  998  999 1000
0  497  499  500  976  979  980  985  986  987  988  989  990  992  993  994  995  996  997  998  999 1000
### pb with 0 if parameter small + 0 always present?

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



