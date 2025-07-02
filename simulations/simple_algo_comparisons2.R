library(dust)

n <- 4
data <- dataGenerator_1D(n, parameters = 0, sdNoise = 1, type = "gauss")
pen <- 2*log(n)
res_dust <- dust.1D(data, penalty = pen/2, model = "gauss", method = "detIndex_Eval1")
res_dust

####

# Définition des vecteurs
#v1 <- res_dust$costQ + rev(cumsum(rev(data^2))) + pen
v1 <- c(-pen/2,res_dust$cost[-n]) + pen
v2 <- (n:1)   # toujours positifs
v3 <-rev(cumsum(rev(data)))

# Domaine d'évaluation
x_vals <- seq(min(v3) - 3, max(v3) + 3, length.out = 1000)

# Calcul des q_i(x) pour chaque i
#q_vals <- sapply(1:length(v1), function(i) v1[i] + v2[i] * (x_vals - v3[i])^2)
q_vals <- sapply(1:length(v1), function(i) v1[i] + v2[i] * x_vals^2 - 2*v3[i] * x_vals)


# Calcul de la fonction enveloppe inférieure : min_i q_i(x)
q_min <- apply(q_vals, 1, min)

# Début du tracé
matplot(x_vals, q_vals, type = "l", lty = 1, col = 1:length(v1), lwd = 1.5,
        ylab = expression(q[t]^i * "(" * theta * ")"),
        xlab = expression(theta),
        main = expression(Q[t]*"("*theta*")" == min[i]*q[t]^i*"("*theta*")"),
        ylim = c(min(q_vals), min(q_vals) + 5))


# Pour chaque fonction, on met en gras les zones où elle est minimale
for (i in 1:ncol(q_vals))
{
  is_min <- abs(q_vals[, i] - q_min) < 1e-8  # tolérance numérique
  # On découpe en segments continus
  rle_min <- rle(is_min)
  idx <- cumsum(rle_min$lengths)
  start <- c(1, head(idx + 1, -1))
  end <- idx

  for (j in which(rle_min$values)) {
    seg_x <- x_vals[start[j]:end[j]]
    seg_y <- q_vals[start[j]:end[j], i]
    lines(seg_x, seg_y, col = i, lwd = 5)
  }
}
legend("topleft",
       legend = lapply(0:(length(v1)-1), function(i) bquote(q[t]^.(i))),
       col = 1:length(v1), lty = 1, lwd = 2)

abline(h = min(q_min) + pen/2, lty = 2, lwd = 2, col = "gray40")

res_dust

# On identifie les indices où la colonne 1 est > colonne 3
dim(q_vals)
indices <- which(q_vals[,1] > q_vals[,4])
point <- indices[which.min(q_vals[indices,4])]
indices

points(x_vals[point], q_vals[point,4], col = "red", pch = 19, cex = 1.5)




res_dust
