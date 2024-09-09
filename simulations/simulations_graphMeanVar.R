
library(dust)

n <- 10^3
beta <- 2*log(n)
y <-  dataGenerator_meanVar(chpts = c(n))
y <- c(rep(0, n), rep(10, n))
plot(y)
r1 <- dust.partitioner.meanVar(method = "randIndex_Eval0")$quick(data = y, penalty = beta)
r2 <- dust.partitioner.meanVar(method = "randIndex_Eval4")$quick(data = y, penalty = beta)

r1$changepoints
r1$nb[n]
r1$nb[n]/n * 100
n*(n+1)/2 / sum(r1$nb)
plot(r1$nb)

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

oneSimu <- function(n)
{
  beta <- 2*log(n)
  #y <-  dataGenerator_meanVar(chpts = c(round(1/3*n), round(2/3*n),n),
  #                            means = c(0,0,1),
  #                            sds = c(1, 2, 2))
  y <-  dataGenerator_meanVar(chpts = n)
  res <- dust.partitioner.meanVar(method = "randIndex_Eval0")$quick(data = y, penalty = beta)
  return(res$nb)
}

# Load required libraries
library(ggplot2)

# Parameters
n_steps <- 1000  # Number of time steps
n_simulations <- 500  # Number of random walks

# Simulate 100 random walks with 1000 time steps
dust_simus <- replicate(n_simulations, oneSimu(n_steps))

beta <- 2*log(n_steps)
y <-  dataGenerator_meanVar(chpts = n_steps)
res <- dust.partitioner.meanVar(method = "randIndex_Eval3")$quick(data = y, penalty = beta)


# Calculate the mean and confidence intervals for each time step
means <- apply(dust_simus, 1, mean)
std_devs <- apply(dust_simus, 1, sd)
up <- apply(dust_simus, 1, function(x) quantile(x,0.95))
down <- apply(dust_simus, 1, function(x) quantile(x,0.05))

# 95% Confidence intervals (mean ± 1.96 * standard error)
#upper_bound <- means + 1.96 * std_devs / sqrt(n_simulations)
#lower_bound <- means - 1.96 * std_devs / sqrt(n_simulations)

upper_bound <- up
lower_bound <- down

# Create a data frame for plotting
time <- 1:n_steps
df <- data.frame(time = time,
                 exemplarOP = c(res$nb),  # Pick the first simulation as the exemplar
                 exemplar = dust_simus[,2],  # Pick the first simulation as the exemplar
                 mean = means,
                 lower = lower_bound,
                 upper = upper_bound)

# Plot using ggplot2
ggplot(df, aes(x = time)) +
  geom_line(aes(y = exemplarOP), color = "green", size = 0.5) +  # Exemplar random walk
  geom_line(aes(y = exemplar), color = "blue", size = 0.5) +  # Exemplar random walk
  geom_line(aes(y = mean), color = "red", size = 1) +  # Mean curve
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +  # Confidence interval
  labs(title = "Random Walk Simulation with Mean and Confidence Interval",
       x = "Time Steps",
       y = "Position",
       caption = "Exemplar random walk (blue, dashed)\nMean curve (red) with 95% confidence interval (grey)") +
  theme_minimal()


############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

n <- 10^4
beta <- 2*log(n)
y <-   dataGenerator_meanVar(chpts = c(n))
plot(y)
system.time(dust.partitioner.meanVar(method = "randIndex_Eval0")$quick(data = y, penalty = beta))
system.time(dust.partitioner.meanVar(method = "randIndex_Eval2")$quick(data = y, penalty = beta))
#system.time(dust.partitioner.meanVar(method = "randIndex_Eval3")$quick(data = y, penalty = beta))
### OP = randIndex_Eval3






