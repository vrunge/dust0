
#library(dust)
#library(changepoint)

#n <- 5*10^4
#y <- rnorm(n)
#system.time(cpt.mean(y, penalty = "Manual", pen.value = 2 * log(n), method = "PELT"))
#system.time(dust.partitioner.1D()$quick(data = y, penalty = log(n)))
#system.time(dust.partitioner.1D(method = "detIndex_Eval5")$quick(data = y, penalty = log(n)))


#plot(ansmean,cpt.col='blue')

#print(ansmean)

# change in mean and variance
#z=c(rnorm(100,0,1),rnorm(100,2,10))
#ansmeanvar=cpt.meanvar(z)
#plot(ansmeanvar,cpt.width=3)
#print(ansmeanvar)
