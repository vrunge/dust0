

#dataGenerator

dataGenerator(chpts = c(50,100), parameter = c(0,1), type = "gauss")
plot(dataGenerator(chpts = c(50,100), parameter = c(0,1), type = "gauss"))


plot(dataGenerator(chpts = c(30,100,120), parameter = c(1,10,1), type = "poisson"))
data <- dataGenerator(chpts = c(300,1000,1200), parameter = c(1,1.4,1), type = "gauss")

\code{"gauss"}, \code{"exp"}, \code{"poisson"}, \code{"geom"}, \code{"bern"}, \code{"binom"}, \code{"negbin"}

type <- "negbin"
s1 <- 12
s2 <- 8
t <- 15
x <- dataGenerator_1D(chpts = 30, parameters = 0.4, type = type)
x <- x/max(x)
mu = seq(0,1,length.out = 1000)
res <- dual_1D(mu = mu, x = x, s1 = s1, s2 = s2, t = t, type = type)
plot(mu, res, type = 'l')

res

(t-s1)/(s1-s2)
