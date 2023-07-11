

#dataGenerator

plot(dataGenerator(chpts = c(30,100,120), parameter = c(1,10,1), type = "gauss"))
plot(dataGenerator(chpts = c(30,100,120), parameter = c(1,10,1), type = "poisson"))



data <- dataGenerator(chpts = c(300,1000,1200), parameter = c(1,1.4,1), type = "gauss")

dust_R(data, 2*log(120))
