library(dust)

n <- 5
data <- dataGenerator_1D(c(n), parameters = c(0), sdNoise = 1, type = "gauss")
#data <- 1:n
#data <- c(-0.5044359, -0.8508462, -1.4807377, 1.0516349,
#          -0.4346470, 2.0146249, -0.5392933, 0.8239567)
pen <- 2*log(n)
res_dust <- dust.1D(data, penalty = pen, model = "gauss", method = "detIndex_Eval1")

level <-  0
####

  # Définition des vecteurs
  #v1 <- res_dust$costQ + rev(cumsum(rev(data^2))) + pen
  v1 <- c(-pen, res_dust$costQ)[1:n] + sum(data[1:n]^2)/2 + pen
  v2 <- (n:1)   # toujours positifs
  v3 <-rev(cumsum(rev(data[1:n])))

  # Domaine d'évaluation
  x_vals <- seq(-2, 2, length.out = 1000)

  # Calcul des q_i(x) pour chaque i
  #q_vals <- sapply(1:length(v1), function(i) v1[i] + v2[i] * (x_vals - v3[i])^2)
  q_vals <- sapply(1:length(v1), function(i) v1[i] + v2[i] * x_vals^2/2 - v3[i] * x_vals)


  # Calcul de la fonction enveloppe inférieure : min_i q_i(x)
  q_min <- apply(q_vals, 1, min)

  bquote(q[.(k)]^.(i))

  # Début du tracé
  matplot(x_vals, q_vals, type = "l", lty = 1, col = 1:length(v1), lwd = 1.5,
          ylab = NA,
          xlab = NA,
          main = bquote(Q[.(n)]*"("*theta*")" == min[i ==0*","*"..."*","*.(n - 1)]*q[.(n)]^i*"("*theta*")"),
          xlim = c(-2, 2),
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
  legend("topleft",
         legend = lapply(0:(length(v1)-1), function(i) bquote(q[.(n)]^.(i))),
         col = 1:length(v1), lty = 1, lwd = 2)

  #abline(h = level + pen, lty = 2, lwd = 2, col = "gray40")
  level <- min(q_min)
  abline(h = level, lty = 3, lwd = 2, col = "gray40")

  # On identifie les indices où la colonne 1 est > colonne 3
  dim(q_vals)
  indices <- which(q_vals[,1] > q_vals[,4])
  point <- indices[which.min(q_vals[indices,4])]
  indices

  points(x_vals[point], q_vals[point,4], col = "red", pch = 19, cex = 1.5)


#plot(q_min)

res_dust
