

### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###

data_gauss_noPruning <- function(n, beta)
{
  sqrt(beta/n) * (sqrt(n-1) - sqrt((1:n)*(n-(1:n))) + sqrt(((1:n)-1)*(n-(1:n)+1)))
}

### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###


data_poisson_noPruning <- function(n, beta, Y)
{
  poisson_zero <- function(t, n, x, Y, beta)
  {(t/n) * x * log(x)  - Y*log(Y) + (Y - (t/n)*x)*log((Y-(t/n)*x)/(1-(t/n))) - beta/n}

  y <- NULL
  for(t in 1:(n-1))
  {
    poisson_zero_x <- function(x) poisson_zero(t, n, x, Y, beta)
    res <- uniroot(poisson_zero_x, interval = c(0.0001, Y))
    y <- c(y, res$root)
  }
  y <- c(y, Y)
  res <- diff((1:n) * y)
  res <- c(y[1], res)

  return(res)
}


### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###


data_exponential_noPruning <- function(n, beta, Y)
{
  exponential_zero <- function(t, n, x, Y, beta)
  {- log((Y-(t/n)*x)/(Y-(t/n)*Y)) + (t/n)*log((Y-(t/n)*x)/(x-(t/n)*x)) - beta/n}

  y <- NULL
  for(t in 1:(n-1))
  {
    exponential_zero_x <- function(x) exponential_zero(t, n, x, Y, beta)
    res <- uniroot(exponential_zero_x, interval = c(0, Y))
    y <- c(y, res$root)
  }
  y <- c(y, Y)
  res <- diff((1:n) * y)
  res <- c(y[1], res)

  return(res)
}



### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###


n <- 30
Y <- 1
### to get mean = Y
beta <- n/(n-1)*Y^2

d1 <- data_gauss_noPruning(n, beta = beta)
d2 <- data_exponential_noPruning(n, beta = 0.1, Y = Y)
d3 <- data_exponential_noPruning(n, beta = 1, Y = Y)


ylimi <- c(min(d1,d2,d3), max(d1,d2,d3))
plot(d1, type = 'b', lwd = 3, ylim = ylimi, cex.axis = 1.5)
par(new = TRUE)
#plot(data_poisson_noPruning(50, 2, 3), type = 'b', lwd = 3, col = 3)
#par(new = TRUE)
plot(d2, type = 'b', lwd = 3, col = 3, ylim = ylimi, cex.axis = 1.5)
par(new = TRUE)
plot(d3, type = 'b', lwd = 3, col = 2,  ylim = ylimi, cex.axis = 1.5)
legend("topleft", legend = c("Gauss", "exponential, beta = 0.1",  "exponential, beta = 1"), col = c(1,3,2), lty = 1, lwd = 3, cex = 2)
mean(d1)
mean(d2)
mean(d3)



### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###

### verification

### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###

resG <- data_gauss_noPruning(50, 2)
resP <- data_poisson_noPruning(50, 2, 3)
resE <- data_exponential_noPruning(50, 2, 1)

plot(resG)
cumy <- ceiling(cumsum(resPoisson))


plot_data_poisson_noPruning <- function(res)
{
  n <- length(res)
  cumy <- cumsum(res)
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

    curve(Cost, from = log(cumy[n]/n)-0.1, to = log(cumy[n] - cumy[n-1])+0.3, ylim = c(mini-0.1,mini+2), n = 1000, col = t)
    par(new = TRUE)
  }

  par(new = FALSE)
}


plot_data_poisson_noPruning(resP)


### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###




plot_data_exponential_noPruning <- function(res)
{
  n <- length(res)
  cumy <- cumsum(res)
  u <- cumy[1]
  res <- 1 + log(u)
  x <- 1/((cumy[n] - cumy[1])/(n-1))
  mini <- (n-1)*(-log(x)) + (cumy[n] - cumy[1])*x + res


  for(t in 1:(n-1))
  {
    Cost <- function(x)
    {
      u <- cumy[t]/t
      res <- t*(1 + 1 *log(u))
      return(-(n-t)*log(x) + (cumy[n] - cumy[t])*x + res)
    }

    curve(Cost, n = 1000, col = t, from = 0, to = 1, ylim = c(mini-0.01, mini+3))
    par(new = TRUE)
  }
  par(new = FALSE)
}


plot_data_exponential_noPruning(resE)



