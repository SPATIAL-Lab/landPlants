
co2Xv = fm(CO2_a = rep(seq(250, 1500, length = 25), 4), 
   v = c(rep(50 * 1e3^2, 25), rep(100 * 1e3^2, 25), rep(200 * 1e3^2, 25), 
         rep(400 * 1e3^2, 25)))

co2Xv$CO2 = rep(seq(250, 1500, length = 25), 4)
co2Xv$v = c(rep(50 * 1e3^2, 25), rep(100 * 1e3^2, 25), rep(200 * 1e3^2, 25), 
            rep(400 * 1e3^2, 25))

layout(matrix(1:4, nrow = 2))
par(mar = c(4, 4, 1, 1))
cats = unique(co2Xv[, 7])
for(i in 1:4){
  dsub = co2Xv[co2Xv[, 7] == cats[1],]
  plot(dsub[, 6], dsub[, 1 + i], xlab = names(dsub)[6], 
       ylab = names(dsub)[1 + i], type = "l", ylim = range(co2Xv[, 1 + i]))
  for(j in 2:length(cats)){
    dsub = co2Xv[co2Xv[, 7] == cats[j],]
    lines(dsub[, 6], dsub[, 1 + i], xlab = names(dsub)[6], 
          ylab = names(dsub)[1 + i], col = j)
  }
  legend("topleft", legend = paste("v =", cats), col = seq_along(cats), 
         lty = 1, bty = "n", bg = NULL)
}
