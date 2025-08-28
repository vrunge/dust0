###
###
### Simon QUERNE. September 2024
###
###

library(dust0)

#######################
### 1. Pruning capacity
#######################

## different cost
## different beta
## different max dual methods
## different signals
## different data size

oneSimu_nb <- function(n, type = "poisson", param = 10, method = "randIndex_Eval0")
{
  beta <- 2*log(n)
  y <-  dataGenerator_1D(chpts = n, parameters = param, type = type)
  res <- dust.partitioner.1D(method = method, model = type)$quick(data = y, penalty = beta)
  cat("0")
  return(res$nb)
}


#####################################################################
library(ggplot2)
# Parameters
n <- 10^5  # Number of time steps
n_simulations <- 1000  # Number of random walks


dust_simus <- replicate(n_simulations, oneSimu_nb(n, method = "fastest"))

#dim(dust_simus)
#df_simus <- as.data.frame(dust_simus)
#dim(df_simus)
#write.csv(df_simus, "your_dataframe.csv", row.names = FALSE)



# Calculate the mean and confidence intervals for each time step
means <- apply(dust_simus, 1, mean)
up <- apply(dust_simus, 1, function(x) quantile(x,0.975))
down <- apply(dust_simus, 1, function(x) quantile(x,0.025))

# Create a data frame for plotting
time <- 1:n
df <- data.frame(time = time,
                 examplar = dust_simus[,1],  # Pick the first simulation as the exemplar
                 mean = means,
                 lower = down,
                 upper = up)
y_min <- min(df$lower)  # You can set a custom value if needed, e.g., -20
y_max <- max(df$upper) # You can set a custom value if needed, e.g., 20


ggplot() +
  geom_line(data = df, aes(x = time, y = examplar, color = "one dust"), size = 0.2) +  # Exemplar random walk
  geom_line(data = df, aes(x = time, y = mean, color = "1000 dust mean"), size = 0.5) +  # Mean curve
  geom_ribbon(data = df,aes(x = time, ymin = lower, ymax = upper), fill = "grey60", alpha = 0.4) +  # Confidence interval
  labs(title = "",
       x = "",
       y = "number of indices",
       caption = "",
       color = "") + ylim(y_min, y_max) +
  scale_color_manual(values = c("PELT" = "orange", "one dust" = "blue", "1000 dust mean" = "red"),  breaks = c("PELT", "one dust", "1000 dust mean")) + # Custom colors for the legend
  theme(
    axis.text = element_text(size = 16),            # Increase size of axis numbers
    axis.title.y = element_text(size = 20),         # Increase size of the x-axis label
    legend.text = element_text(size = 20),          # Increase size of the legend text
    legend.title = element_text(size = 16),          # Increase size of the legend title
    legend.position = "top"
  )



# Step 2: Apply log transformation
log_x <- log((1:n)[(n/2):n])
log_y <- log(means[(n/2):n])

# Step 3: Fit a linear model for the log-log data
model <- lm(log_y ~ log_x)

# Step 4: Plot log-log plot with regression line
plot(log_x, log_y, main="Log-Log Plot", xlab="log(x)", ylab="log(y)", pch=19)
abline(model, col="blue")

summary(model)


######################
### 2. Time comparison
######################

## with fpop and gfpop with the fastest max dual method





