


data <- dataGenerator_1D(chpts = 1, parameters = c(0,1), sdNoise = 0.2, type = "gauss")
data <- data_normalization_1D(1)
plot(data)
res <- dust.1D(data = 1, penalty = 30)


#################################################################################
#################################################################################
### TEST dust.1D

data <- dataGenerator_1D(chpts = c(50,100), parameters = c(0,1), sdNoise = 0.2, type = "gauss")
data <- data_normalization_1D(data)
plot(data)
res <- dust.1D(data = data)
res$changepoints
res$lastIndexSet
res$nb
res$costQ

res <- dust.1D(data = data, penalty = 0)
res$changepoints

# 0: random eval
# 1: exact eval (if possible, otherwise, -inf (OP))
# 2: golden-section search
# 3: binary search. At each step, we evaluate the tangent line to the current point at its max to stop the search at early step (when possible)
# 4: Quasi-Newton
# 5: PELT
# 6: OP


data <- dataGenerator_1D(chpts = c(500,600), parameters = c(0,1), sdNoise = 0.2, type = "gauss")
plot(data)
res <- dust.1D(data = data)
COST <- res$costQ


for(i in 0:6)
{
res <- dust.1D(data = data, method = paste0("randIndex_Eval",i))
print(c("randIndex_Eval", i))
print(res$lastIndexSet)
print(res$changepoints)
print(all(res$costQ == COST))
}

for(i in 0:6)
{
  res <- dust.1D(data = data, method = paste0("detIndex_Eval",i))
  print(c("detIndex_Eval", i))
  print(res$lastIndexSet)
  print(all(res$costQ == COST))
}


####
model <- c("gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance")

for(i in 1:7)
{

  data <- dataGenerator_1D(chpts = c(500,1000), parameters = c(0.5,0.9), type = model[i])
  res <- dust.1D(data = data, model = model[i])
  COST <- res$costQ
  for(j in 0:6)
  {
    res <- dust.1D(data = data,
                   model = model[i],
                   method = paste0("randIndex_Eval",j))
    print(c("model ", model[i], " randIndex_Eval ", j))
    print(res$lastIndexSet)
    print(res$changepoints)
    print(all(res$costQ == COST))
  }

}



#################################################################################
#################################################################################
### TEST dust.object.1D


data = data
penalty = 2*log(length(data))
model = "gauss"
method = "fastest"
nbLoops = 10


partitioner <- new(DUST_1D, model, method, nbLoops)
partitioner$one_dust(data, penalty)
partitioner$get_partition()

partitioner <- dust.object.1D()
partitioner$
partitioner$get_partition()





