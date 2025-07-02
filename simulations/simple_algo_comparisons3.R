library(dust)

n <- 7
data <- dataGenerator_1D(c(n), parameters = c(0), sdNoise = 1, type = "gauss")
#data <- 1:n
#data <- c(-0.5044359, -0.8508462, -1.4807377, 1.0516349,
#          -0.4346470, 2.0146249, -0.5392933, 0.8239567)
pen <- 2*log(n)
res_dust <- dust.1D(data, penalty = pen, model = "gauss", method = "detIndex_Eval1")

layout(matrix(1:2, ncol = 1), heights = c(4, 1))

level <-  0
####
for(k in 1:n)
{
  # Définition des vecteurs
  #v1 <- res_dust$costQ + rev(cumsum(rev(data^2))) + pen
  v1 <- c(-pen, res_dust$costQ)[1:k] + sum(data[1:k]^2)/2 + pen
  v2 <- (k:1)   # toujours positifs
  v3 <-rev(cumsum(rev(data[1:k])))

  # Domaine d'évaluation
  x_vals <- seq(-2, 4, length.out = 1000)

  # Calcul des q_i(x) pour chaque i
  #q_vals <- sapply(1:length(v1), function(i) v1[i] + v2[i] * (x_vals - v3[i])^2)
  q_vals <- sapply(1:length(v1), function(i) v1[i] + v2[i] * x_vals^2/2 - v3[i] * x_vals)


  # Calcul de la fonction enveloppe inférieure : min_i q_i(x)
  q_min <- apply(q_vals, 1, min)

  bquote(q[.(k)]^.(i))

  # Début du tracé
  par(mar = c(1, 1, 3, 1))  # marges normales
  matplot(x_vals, q_vals, type = "l", lty = 1, col = 1:length(v1), lwd = 1.5,
          ylab = NA,
          xlab = NA,
          main = bquote(Q[.(k)]*"("*theta*")" == min[i ==0*","*"..."*","*.(k - 1)]*q[.(k)]^i*"("*theta*")"),
          xlim = c(-2, 4),
          ylim = c(- pen/2, 15))


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
  legend("topright",
         legend = lapply(0:(length(v1)-1), function(i) bquote(q[.(k)]^.(i))),
         col = 1:length(v1), lty = 1, lwd = 2)

  #abline(h = level + pen, lty = 2, lwd = 2, col = "gray40")
  level <- min(q_min)
  abline(h = level, lty = 3, lwd = 2, col = "gray40")

  par(mar = c(1, 1, 1, 1))  # marges normales
  eps <- 0.2*(max(data) - min(data))
  plot(x = 1:n, c(data[1:k], rep(NA,n-k)), type = 'b', ylim = c(min(data)-eps,max(data)+eps))
}

#plot(q_min)

res_dust
