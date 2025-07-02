


data <- dataGenerator_meanVar(chpts = c(30,100,120), means = c(0,1,0), sds = c(1,1,2))
res <- dust_R_2param_meanVar(data, 4*log(120), pruningOpt = 2)
res$costQ


# Parameters
test <- -1

while(test < 0)
{
  A <-  runif(1, min = -1, max = 1)
  D <-  runif(1, min = -1, max = 1)
  a <-  runif(1, min = -1, max = 1)
  d <- runif(1, min = -1, max = 1)
  DELTA <- runif(1, min = -1, max = 1)
  A
  D
  a
  d
  DELTA
  test <- A-a^2
}

u-v
u

# Function definition
f <- function(x) {
  inside_log <- A + D * x - (a + d * x)^2
  ifelse(inside_log > 0,
         0.5 * log(inside_log) - DELTA * x,
         NA)  # avoid log of negative or zero
}


u <- (D-2*a*d)/(2*d^2)
v <- sqrt(((D-2*a*d)/(2*d^2))^2 + (A-a^2)/d^2)
u+v
u-v

epsilon <- 2*abs(u+v)/10000

# Plotting range
x_vals <- seq(0, u+v-epsilon, length.out = 1000)
y_vals <- sapply(x_vals, f)

# Plot
plot(x_vals, y_vals, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "f(x)",
     main = "Plot of f(x) = 1/2 * log(A + D x - (a + d x)^2) - DELTA x")
grid()
abline(v = 0, col = 2, lwd =2)

u2 <- (D-2*a*d)/(2*d^2)  + 1/(2*DELTA)
v2 <- sqrt(((D-2*a*d)/(2*d^2))^2 + 1/(4*DELTA^2) + (A-a^2)/d^2)



u2-sign(DELTA)*v2
abline(v = u2-sign(DELTA)*v2, col = 3)
DELTA
u - v < u2-sign(DELTA)*v2
u2-sign(DELTA)*v2 < u + v
u + v
sign(DELTA)
