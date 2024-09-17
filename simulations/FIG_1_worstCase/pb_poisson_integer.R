

n <- 20
Y <- 30


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
}

y <- c(y, Y)
yy <- diff((1:length(y))*(y))
resPoisson <- c(y[1],yy)

mean(resPoisson)


plot(resPoisson, type = 'b',
     ylim = c(min(resPoisson),max(resPoisson)),
     col = 2, lwd = 3)

round(resPoisson)


##################################################################
##################################################################

ceiling(resPoisson)

#cumy <- ceiling(cumsum(resPoisson))
cumy <- round(cumsum(resPoisson))


l <- 0.3
r <- 3
a <- 3.3
b <- 4.5
plot_data_poisson_noPruning <- function(cumy)
{
  u <- cumy[1]
  const <- u - u * log(u)
  if(cumy[1] == 0){const <- 0}
  x <- log((cumy[n] - cumy[1])/(n-1))
  mini <- (n-1)*exp(x) - (cumy[n] - cumy[1])*x + const

  for(t in 1:(n-1))
  {
    Cost <- function(x)
    {
      u <- cumy[t]/t
      const <- t*(u - u *log(u))
      if(cumy[t] == 0){const <- 0}
      return((n-t)*exp(x) - (cumy[n] - cumy[t])*x + const)
    }

    curve(Cost, n = 1000, col = t+1, from = a, to = b, ylim = c(mini-l, mini+r), lwd = 2)
    par(new = TRUE)
  }
  Cost <- function(x)
  {
    return(n*exp(x) - (cumy[n])*x - beta)
  }
  curve(Cost, n = 1000, col = 1, from = a, to = b, ylim = c(mini-l, mini+r), lwd = 2)
  par(new = FALSE)
}


plot_data_poisson_noPruning(cumy)

















