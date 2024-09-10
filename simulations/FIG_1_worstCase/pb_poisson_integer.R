

n <- 40
Y <- 1000
beta <- n/(n-1)*Y^2


sol <- function(t,n, beta)
{
  sqrt(beta/n) * (sqrt(n-1) - sqrt(t*(n-t)) + sqrt((t-1)*(n-t+1)))
}

resGAUSS <- sol(1:n,n,beta)

################################################################################


max <- n/(n-1)*Y
beta <- Y*log(n/(n-1))*9.95/10*n

#### SOLVE = 0
poisson_zero_2 <- function(t, n, x, Y, beta)
{
  (t/n) * x * log(x)  - Y*log(Y) + (Y - (t/n)*x)*log((Y-(t/n)*x)/(1-(t/n))) - beta/n
}

y <- NULL

for(t in 1:(n-1))
{

  poisson_zero_3 <- function(x) poisson_zero_2(t, n, x, Y, beta)

  res <- uniroot(poisson_zero_3, interval = c(0.0001, Y))
  y <- c(y, res$root)
  w <- seq(from = 10^(-10), to = Y, length.out = 100)
  z <- poisson_zero_3(w)
  #plot(w, z, type = 'l', lwd = 1)
  #abline(h = 0, col = 2, lwd = 1)
}

y <- c(y, Y)
yy <- diff((1:length(y))*(y))
resPoisson <- c(y[1],yy)


mean(resGAUSS)
mean(resPoisson)

plot(resGAUSS, type = 'b', ylim = c(min(resGAUSS, resPoisson),max(resGAUSS, resPoisson)), lwd = 3)
par(new = TRUE)
plot(resPoisson, type = 'b', ylim = c(min(resGAUSS, resPoisson),max(resGAUSS, resPoisson)), col = 2, lwd = 3)



round(resPoisson)




cumy <- ceiling(cumsum(resPoisson))

#cumy <- (cumsum(resPoisson))


u <- cumy[1]
res <- (u - u *log(u))
x <- log((cumy[n] - cumy[1])/(n-1))
mini <- (n-1)*exp(x) - (cumy[n] - cumy[1])*x + res




for(t in 1:(n-1))
{
  Cost <- function(x)
  {
    u <- cumy[t]/t
    res <- t*(u - u *log(u))
    return((n-t)*exp(x) - (cumy[n] - cumy[t])*x + res)
  }

  curve(Cost, from = log(cumy[n]/n)-0.1, to = log(cumy[n] - cumy[n-1])+0.1, ylim = c(mini-5,mini+100), n = 1000, col = t)
  par(new = TRUE)
}



