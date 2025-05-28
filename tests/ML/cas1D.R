
getwd()

# nombre de données à générer
n <- 10000
library(dust)

#####
##### GENERATE DATA
#####


###
### gauss
###

epsilon <- 0.3
data <- dataGenerator_1D(chpts = (1:10)*n/10, parameters = 1:10*epsilon)
plot(data)
dust.1D(data, method = "detIndex_Eval6", model = "gauss")
data <- c(rep(0,1000),rep(1,1000))
dust.1D(data, method = "detIndex_Eval6", model = "gauss")

###
### poisson
###

data <- dataGenerator_1D(chpts = n, parameters = 5, type = "poisson")
plot(data)
dust.1D(data, method = "detIndex_Eval6", model = "poisson")
data <- dataGenerator_1D(chpts = n, parameters = 0.5, type = "poisson")
plot(data)
dust.1D(data, method = "detIndex_Eval6", model = "poisson")


###
### exp
###

data <- dataGenerator_1D(chpts = n, parameters = 50, type = "exp")
plot(data)
dust.1D(data, method = "detIndex_Eval6", model = "exp")
data <- dataGenerator_1D(chpts = n, parameters = 0.01, type = "exp")
plot(data)
dust.1D(data, method = "detIndex_Eval6", model = "exp")


###
### geom
###

data <- dataGenerator_1D(chpts = n, parameters = 0.01, type = "geom")
plot(data)
res <- dust.1D(data, method = "detIndex_Eval6", model = "geom")
data <- dataGenerator_1D(chpts = n, parameters = 0.95, type = "geom")
plot(data)
res <- dust.1D(data, method = "detIndex_Eval6", model = "geom")



###
### bern
###

data <- dataGenerator_1D(chpts = n, parameters = 0.9, type = "bern")
plot(data)
res <- dust.1D(data, method = "detIndex_Eval6", model = "bern")
data <- dataGenerator_1D(chpts = n, parameters = 0.05, type = "bern")
plot(data)
res <- dust.1D(data, method = "detIndex_Eval6", model = "bern")

###
### binom
###

data <- dataGenerator_1D(chpts = n, parameters = 0.95, type = "binom", nbTrials = 100)
plot(data)
data <- data_normalization_1D(data, type = "binom")
res <- dust.1D(data, method = "detIndex_Eval6", model = "binom")

data <- dataGenerator_1D(chpts = n, parameters = 0.05, type = "binom")
data <- data_normalization_1D(data, type = "binom")
res <- dust.1D(data, method = "detIndex_Eval6", model = "binom")



###
### negbin
###

data <- dataGenerator_1D(chpts = n, parameters = 0.3, type = "negbin", nbSuccess = 30)
plot(data)
data <- data_normalization_1D(data, type = "negbin")
plot(data)
res <- dust.1D(data, method = "detIndex_Eval6", model = "negbin")

data <- dataGenerator_1D(chpts = n, parameters = 0.95, type = "negbin")
plot(data)
data <- data_normalization_1D(data, type = "negbin")
plot(data)
res <- dust.1D(data, method = "detIndex_Eval6", model = "negbin")



###
### variance
###

epsilon <- 0.3
data <- dataGenerator_1D(chpts = (1:10)*n/10, parameters = 1:10*epsilon, type = "variance")
plot(data)
dust.1D(data, method = "detIndex_Eval6", model = "variance")
data <- c(rep(0,1000),rep(1,1000))
dust.1D(data, method = "detIndex_Eval6", model = "variance")

#####
##### LOAD DATA
#####

data_gauss <- read_csv("dataset_1D_gauss.csv")
data_poisson <- read_csv("dataset_1D_poisson.csv")
data_exp <- read_csv("dataset_1D_exp.csv")
data_geom <- read_csv("dataset_1D_geom.csv")
data_bern <- read_csv("dataset_1D_bern.csv")
data_binom <- read_csv("dataset_1D_binom.csv")
data_negbin <- read_csv("dataset_1D_negbin.csv")
data_variance <- read_csv("dataset_1D_variance.csv")



#####
##### DATA ANALYSIS
#####

sum(is.na(data_gauss))
sum(is.na(data_poisson))
sum(is.na(data_exp))
sum(is.na(data_geom))
sum(is.na(data_bern))
sum(is.na(data_binom))
sum(is.na(data_negbin))
sum(is.na(data_variance))

table(data_gauss$ab)
table(data_poisson$ab)
table(data_exp$ab)
table(data_geom$ab)
table(data_bern$ab)
table(data_binom$ab)
table(data_negbin$ab)
table(data_variance$ab)

c1 <- "a!=b"
c2 <- "a=b"
c3 <- "f-nf"
c4 <- "f-f"

table(data_gauss$pruning[data_gauss$ab == c1])
table(data_gauss$pruning[data_gauss$ab == c2])


table(data_poisson$pruning[data_poisson$ab == c1])
table(data_poisson$pruning[data_poisson$ab == c2])
table(data_poisson$pruning[data_poisson$ab == c3])
table(data_poisson$pruning[data_poisson$ab == c4])

table(data_geom$pruning[data_geom$ab == c1])
table(data_geom$pruning[data_geom$ab == c2])
table(data_geom$pruning[data_geom$ab == c3])
table(data_geom$pruning[data_geom$ab == c4])

table(data_bern$pruning[data_bern$ab == c1])
table(data_bern$pruning[data_bern$ab == c2])
table(data_bern$pruning[data_bern$ab == c3])
table(data_bern$pruning[data_bern$ab == c4])

table(data_binom$pruning[data_binom$ab == c1])
table(data_binom$pruning[data_binom$ab == c2])
table(data_binom$pruning[data_binom$ab == c3])
table(data_binom$pruning[data_binom$ab == c4])

table(data_negbin$pruning[data_negbin$ab == c1])
table(data_negbin$pruning[data_negbin$ab == c2])
table(data_negbin$pruning[data_negbin$ab == c3])
table(data_negbin$pruning[data_negbin$ab == c4])




data_gauss2 <- data_gauss[data_gauss$ab == c1,]
data_poisson2 <- data_poisson[data_poisson$ab == c1,]
data_exp2 <-  data_exp[data_exp$ab == c1,]
data_geom2 <- data_geom[data_geom$ab == c1,]
data_bern2 <- data_bern[data_bern$ab == c1,]
data_binom2 <- data_binom[data_binom$ab == c1,]
data_negbin2 <- data_negbin[data_negbin$ab == c1,]
data_variance2 <- data_variance[data_variance$ab == c1,]


hist(data_gauss2$mu)

data_poisson2$muMax - pmin(1,(data_poisson2$objectiveMean)/(data_poisson2$constraintMean))

plot(data_geom2$objectiveMean/data_geom2$constraintMean)

data_gauss2 <- cbind(data_gauss2,data_gauss2$constraintMean/data_gauss2$objectiveMean)
data_gauss2 <- cbind(data_gauss2,data_gauss2$/data_gauss2$objectiveMean)


data_gauss <- read_csv("dataset_1D_gauss.csv")
write_csv(data_gauss2, file = "dataset_1D_gauss2.csv")

