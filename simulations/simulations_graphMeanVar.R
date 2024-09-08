

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
# Load ggplot2 library
library(ggplot2)

# Example data: two vectors of the same length
x <- 1:n # X-axis (e.g., time or indices)
y1 <- r1$nb                        # First vector (e.g., sine wave)
y2 <- r2$nb                    # Second vector (e.g., cosine wave)

# Create a data frame
df <- data.frame(x = x, y1 = y1, y2 = y2)

# Reshape the data to long format for ggplot
df_long <- reshape2::melt(df, id.vars = "x")

# Plot using ggplot2
ggplot(df_long, aes(x = x, y = value, color = variable)) +
  geom_line(size = 1) +
  labs(title = "Plot of Two Vectors",
       x = "",
       y = "number of indices",
       color = "Vectors") +        # Legend title
  theme_minimal()





############################################################################
############################################################################
n <- 10^5
beta <- 2*log(n)
y <-   dataGenerator_meanVar(chpts = c(n))
plot(y)
system.time(dust.partitioner.meanVar(method = "randIndex_Eval0")$quick(data = y, penalty = beta))
system.time(dust.partitioner.meanVar(method = "randIndex_Eval4")$quick(data = y, penalty = beta))
system.time(dust.partitioner.meanVar(method = "randIndex_Eval3")$quick(data = y, penalty = beta))
### OP = randIndex_Eval3

