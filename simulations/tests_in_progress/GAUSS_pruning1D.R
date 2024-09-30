

t <- 10^3

y <- rnorm(t)
cumy <- c(0,cumsum(y))

r <- 20
s <- 50

meany <- function(u,v)
{
 (cumy[v+1] - cumy[u+1])/(v-u)
}
Qdiff <- function(u,v)
{
  (v*meany(0,v)^2 - u*meany(0,u)^2) / (v - u)
}

A <- -meany(r,s)^2 + Qdiff(r,s)
B <- meany(r,s)*meany(s,t) - Qdiff(r,s) - Qdiff(s,t)
C <- -meany(s,t)^2 + Qdiff(s,t)


curve(A*x^2 + B*x+C, xlim = c(0,1), ylim = c(-0.1,0.1))
abline(v=0,h=0)


###########################

z <- NULL
m <- 10^3
n <- 200
for(i in 1:m)
{
  y <- rnorm(n)
  z <- c(z, sqrt(sum(y^2)))
}

hist(z, breaks = 50, probability = TRUE,
     main = "Histogram with Bell Curve",
     xlab = "z values", col = "lightblue")

# Add the bell curve (normal distribution)
curve(dnorm(x, mean = base::mean(z), sd = sd(z)),
      col = "red", lwd = 2, add = TRUE)



#################################################################################
#################################################################################


t <- 101

y <- rnorm(t)
cumy <- c(0,cumsum(y))

s <- 100
r <- 1


################################################################################################
################################################################################################
################################################################################################

pruning <- function(r,s,t)
{
  y <- rnorm(t)
  cumy <- c(0,cumsum(y))
  meany <- function(u,v){(cumy[v+1] - cumy[u+1])/(v-u)}
  res <- sqrt(r/s)*abs(meany(0,r) - meany(r,s))  - abs(meany(s,t) - meany(r,s)) - sqrt(s/t)*abs(meany(0,s) - meany(s,t))
  return(res)
}

################################################

pruning_dim <- function(r,s,t, dim = 2)
{
  y <- matrix(rnorm(t * dim), nrow = dim, ncol = t)
  cumy <- cbind(0,t(apply(y, 1, cumsum)))
  meany <- function(u,v){(cumy[,v+1] - cumy[,u+1])/(v-u)}
  res <- sqrt(r/s)*sqrt(sum((meany(0,r) - meany(r,s))^2)) - sqrt(sum((meany(s,t) - meany(r,s))^2)) - sqrt(s/t)*sqrt(sum((meany(0,s) - meany(s,t))^2))
  return(res)
}

################################################

m <- 10^5
t <- 10
s <- 9
r <- 8

res <- replicate(m, pruning(r,s,t))
hist(res, breaks = 200)
abline(v = 0, col = 2, lwd =2)
mean(res>0)


res <- replicate(m, pruning_dim(r,s,t, dim = 4))
hist(res, breaks = 200, probability = TRUE)
abline(v = 0, col = 2, lwd =2)
mean(res>0)

curve(dnorm(x, mean = mean(res), sd = sd(res)),
      col = "blue", lwd = 2, add = TRUE)




val <- NULL
for(i in 1:10)
{
  print(i)
  res <- replicate(m, pruning_dim(r,s,t, dim = i))
  val <- c(val, mean(res>0))
}

val
plot(val)
plot((-log(val)))
plot((1:10), (-log(val)))



################################################################################################
################################################################################################
################################################################################################

pruningALL <- function(s, t)
{
  y <- rnorm(t)
  cumy <- c(0,cumsum(y))
  meany <- function(u,v){(cumy[v+1] - cumy[u+1])/(v-u)}
  res <- rep(0, s-1)
  for(r in 1:(s-1))
  {
    res[r] <- sqrt(r/s)*abs(meany(0,r) - meany(r,s))  - abs(meany(s,t) - meany(r,s)) - sqrt(s/t)*abs(meany(0,s) - meany(s,t))
  }
  return(res)
}

res <- pruningALL(30,100)
plot(res, type = 'b', ylim =c(min(c(res,-0.01)), max(c(res, 0.01))))
abline(h=0, col = 2, lwd = 2)

################################################



m <- 10^5
t <- 10000
s <- 50
r <- 8

res <- replicate(m, pruning(r,s,t))
hist(res, breaks = 200, probability = T)
abline(v = 0, col = 2, lwd =2)
mean(res>0)

curve(dnorm(x, mean = mean(res), sd = sd(res)),
      col = "blue", lwd = 2, add = TRUE)


#################################################


pruningALLALL <- function(t)
{
  y <- rnorm(t)
  cumy <- c(0,cumsum(y))
  meany <- function(u,v){(cumy[v+1] - cumy[u+1])/(v-u)}

  res <- matrix(rep(-Inf, (t-1)*(t-1)), t-1, t-1)
  for(s in 2:(t-1))
  {
    for(r in 1:(s-1))
    {
      res[r,s] <- sqrt(r/s)*abs(meany(0,r) - meany(r,s))  - abs(meany(s,t) - meany(r,s)) - sqrt(s/t)*abs(meany(0,s) - meany(s,t))
    }
  }
  return(res)
}

t <- 1000
res <- pruningALLALL(t)
#res[t-1,1:(t-1)/2] <- 1 #positif marron
res2 <- t(res[nrow(res):1,])

image(res2)
image(res2>=0)

apply(res, 2, max)
t - sum(apply(res, 2, function(x) max(x) > 0))


res1 <- dust.partitioner.1D(model = "gauss", method = "fastest")$quick(data = rnorm(t), penalty = log(t))
res1$changepoints
res1$nb





########################################################################
########################################################################






PruningIndex_s_time_t <- function(s, t)
{
  y <- rnorm(t)
  cumy <- c(0,cumsum(y))
  meany <- function(u,v){(cumy[v+1] - cumy[u+1])/(v-u)}

  res <- matrix(rep(-Inf, (s-1)*(t-s)), s-1, t-s)
  for(u in 1:(t-s))
  {
    tt <- s + u
    for(r in 1:(s-1))
    {
      res[r,u] <- sqrt(r/s)*abs(meany(0,r) - meany(r,s))  - abs(meany(s,tt) - meany(r,s)) - sqrt(s/tt)*abs(meany(0,s) - meany(s,tt))
    }
  }
  return(res)
}

s <- 50
t <- 200
res <- PruningIndex_s_time_t(s, t)
res
#res[t-1,1:(t-1)/2] <- 1 #positif marron
res2 <- t(res[nrow(res):1,])

image(res2)
image(res2>=0)

which(apply(res, 2, max) > 0)[1]
max(apply(res, 2, max)) > 0

