

###
### data from dataGenerator_1D
### plots for the readme file
###

par(mfrow=c(1,2))
data <- dataGenerator_1D(chpts = c(50,80, 100), parameters = c(0,1,-1), sdNoise = 0.2, type = "gauss")
plot(data, type = 'b', xlab = "GAUSS COST with 'chpts = c(50,80,100), parameters = c(0,1,-1), sdNoise = 0.2'", ylab = "", col = 4)

#data <- dataGenerator_1D(chpts = c(50,100), parameters = c(10,20), gamma = c(0.9,0.95), type = "gauss")
#plot(data, type = 'b', xlab = "GAUSS COST with 'chpts = c(50,100), parameters = c(10,20), gamma = c(0.9,0.95)'", ylab = "")

data <- dataGenerator_1D(chpts = c(40,100,130), parameters = c(2,14,5), type = "exp")
plot(data, type = 'b', xlab = "EXP COST with 'chpts = c(40,100,130), parameters = c(2,14,5)'", ylab = "", col = 4)


###################################


par(mfrow=c(2,2))
par(mar = c(5, 3, 1, 1))
data <- dataGenerator_1D(chpts = c(50,100), parameters = c(3,10), type = "poisson")
plot(data, type = 'b', xlab = "POISSON COST with 'chpts = c(50,100), parameters = c(3,10)'", ylab = "", col = 4)

data <- dataGenerator_1D(chpts = c(70,120,200), parameters = c(0.7,0.2,0.5), type = "geom")
plot(data, type = 'b', xlab = "GEOM COST with 'chpts = c(70,120,200), parameters = c(0.7,0.3,0.6)'", ylab = "", col = 4)

data <- dataGenerator_1D(chpts = c(30,80,110), parameters = c(0.7, 0.1, 0.3), nbTrials = 5, type = "binom")
plot(data, type = 'b', xlab = "BINOM COST with 'chpts = c(30,80,110), parameters = c(0.7, 0.1, 0.3), nbTrials = 5'", ylab = "", col = 4)

data <- dataGenerator_1D(chpts = c(50,100), parameters = c(0.4,0.7), nbSuccess = 10, type = "negbin")
plot(data, type = 'b', xlab = "NEGBIN COST with 'chpts = c(50,100), parameters = c(0.4,0.7), nbSuccess = 10'", ylab = "", col = 4)





