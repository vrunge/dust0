


mu <- c(rep(0,2500), rep(0.1, 2500))
data <- matrix(rnorm(1000 * 5, mu), nrow = 5, ncol = 1000, byrow = FALSE)
plot(data[3,])

pen <- 2*5*log(1000)/10

res1 <- dust.MD(data, method = "randIndex_Eval4",
               penalty = pen, constraints_l = 1)
res1$changepoints
res1$nb
###

res2 <- dust.MD(data, method = "randIndex_Eval3",
               penalty = pen, constraints_l = 1)
res2$changepoints
res2$nb
##
all(res1$costQ == res2$costQ)

################################################
################################################
################################################

d <- 1
mu <- c(rep(0,250), rep(2, 250))
data <- matrix(rnorm(500* d, mu), nrow = 500, ncol = d, byrow = FALSE)
data <- t(data)
plot(data[1,])

pen <- 2*d*log(500)



res1 <- dust.MD(data, method = "randIndex_Eval1",
                penalty = pen, constraints_l = 1)
res1$changepoints
res1$nb
###


res2 <- dust.MD(data, method = "randIndex_Eval0",
                penalty = pen, constraints_l = 1)
res2$changepoints
res2$nb
res2bis <- dust.1D(as.vector(data), method = "randIndex_Eval0",
                   penalty = pen)
res2bis$changepoints
res2bis$nb
###


d <- 5
mu <- c(rep(0,250), rep(10, 250))
data <- matrix(rnorm(500* d, mu), nrow = 500, ncol = d, byrow = FALSE)
data <- t(data)
plot(data[1,])

pen <- 2*d*log(500)


res3 <- dust.MD(data,
                method = "randIndex_Eval0",
                penalty = pen,
                constraints_l = 1,
                constraints_r = 1)
res3$changepoints
res3$nb
###

res4 <- dust.MD(data,
                method = "randIndex_Eval1",
                penalty = pen,
                constraints_l = 1,
                constraints_r = 2)
res4$changepoints
res4$nb

res4$nb - res3$nb

###

all(res3$costQ == res4$costQ)


#############################################




res4 <- dust.MD(data, method = "randIndex_Eval3",
                penalty = pen, constraints_l = 2, constraints_r = 2)
res4$changepoints
res4$nb



#############################################
#############################################
#############################################

d <- 5
cl <- 1
cr <- 4
n <- 500
mu <- c(rep(4,n/2), rep(4, n/2))
data <- matrix(rnorm(n* d, mu), nrow = n, ncol = d, byrow = FALSE)
data <- t(data)
#plot(data[1,])
dim(data)
pen <- 2*d*log(n)


res0 <- dust.MD(data,
                method = "detIndex_Eval0",
                penalty = pen,
                constraints_l = cl,
                constraints_r = cr)
res0$changepoints
res0$nb
res0$nb[length(res0$nb)]
###

######
###### HERE PB WHEN detIndex_Eval1 AND cr = d-1
######
res1 <- dust.MD(data,
                method = "detIndex_Eval1",
                penalty = pen,
                constraints_l = cl,
                constraints_r = cr)
res1$changepoints
res1$nb
res1$nb[length(res1$nb)]

res1$nb - res0$nb
#plot(res1$nb - res0$nb)
#abline(h = 0, col = 2)

###


res3 <- dust.MD(data,
                method = "randIndex_Eval3",
                penalty = pen,
                constraints_l = cl,
                constraints_r = cr, nbLoops = d)
res3$changepoints
res3$nb
res3$nb[length(res3$nb)]

res3$nb - res0$nb
#plot(res3$nb - res0$nb)
#abline(h = 0, col = 2)


all(res1$costQ == res0$costQ)
all(res3$costQ == res0$costQ)

#plot(res3$nb - res0$nb)



library(ggplot2)

x <- 1:length(res0$nb)


# Create a data frame in long format manually
df <- data.frame(
  x = rep(x, 3),
  y = c(res0$nb, res1$nb, res3$nb),
  group = factor(rep(c("algo0", "algo1", "algo3"), each = length(res0$nb)))
)

# Plot with ggplot2
ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(size = 0.5) +
  geom_point(size = 0) +
  labs(title = "number of indices to save",
       x = "X-axis",
       y = "Values",
       color = "Legend") +
  theme_minimal()

