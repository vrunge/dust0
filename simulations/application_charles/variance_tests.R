

library(dust0)

n <- 10^3
nb <- 10
n/nb
chgpt <- seq(from = 0, to = 1, length.out = nb)[-1]

D <- 4
df <- as.data.frame(sapply(1:D, function(i) sample(nb - 1)))
colnames(df) <- paste0("ts", 1:D)

df

data <- dataGenerator_MD(chpts = floor(chgpt * n),
        parameters = df,
        type = "variance")
dim(data)

# Configuration du layout : 3 lignes, 1 colonne
par(mfrow = c(1,1))

# Tracé des 3 courbes séparément
for (i in 1:3)
{
  #plot(data[i, ], type = "l", col = i, lwd = 2,
  #     xlab = "Index", ylab = paste("Ligne", i),
  #     main = paste("Courbe", i))
}


system.time(
  res <- dust.MD(data,
                 model = "variance",
                 method = "detIndex_Eval5",
                 constraints_l = 1,
                 constraint_r = 0,
                 nbLoops = 10
                )
          )[[1]]

res$changepoints
floor(chgpt * n)
plot(res$nb)




system.time(
  res2 <- dust.MD(data,
                  model = "variance",
                  method = "detIndex_Eval5",
                  nbLoops = 10
  )
)[[1]]

res2$changepoints
floor(chgpt * n)
plot(res2$nb)


res <- dust.1D(dataGenerator_1D(c(500,1000), parameters = c(0,1)), penalty = 2*log(1000), model = "gauss")
plot(res$nb)

