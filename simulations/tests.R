

#dataGenerator

dataGenerator(chpts = c(50,100), parameter = c(0,1), type = "gauss")
plot(dataGenerator(chpts = c(50,100), parameter = c(0,1), type = "gauss"))


plot(dataGenerator(chpts = c(30,100,120), parameter = c(1,10,1), type = "poisson"))
data <- dataGenerator(chpts = c(300,1000,1200), parameter = c(1,1.4,1), type = "gauss")

\code{"gauss"}, \code{"exp"}, \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}

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
