library(R2jags)
library(openxlsx)
source("code/helpers.R")

# To save
parms = c("Pl", "l", "amax.scale", "D", "gc.scale", "ca", "meso.scale",
          "Ci0_m", "A0_m", "d13Ca_m", "A", "D13C", "gcop")

# Zhang data, multi-site, multi-species ----
d = read.xlsx("data/stomata-franks_zhang_2024_P1.0.xlsx", sheet = 1, startRow = 3)

## First w/o collapsing, single specimen interpretations
data = parseFranks(d, FALSE)
system.time({
  post = jags.parallel(data, NULL, parms, "code/models/forwardFranksMulti.R", n.chains = 3,
                       n.iter = 2e6, n.burnin = 1e4, n.thin = 1e3)
})

## Now w/ collapsing
data = parseFranks(d)
system.time({
  post.cl = jags.parallel(data, NULL, parms, "code/models/forwardFranksMulti.R", n.chains = 3,
                          n.iter = 2e6, n.burnin = 1e4, n.thin = 1e3)
})

View(post$BUGSoutput$summary)
View(post.cl$BUGSoutput$summary)

## Plot comparing sample medians for PSM and trad
plot(d$CO2_ppm[-(1:3)], post$BUGSoutput$median$ca, pch = 21, bg = "grey50",
     cex = 2, lwd = 2, ylab = "Median ca from PSM", 
     xlab = "Median ca traditional")
abline(0, 1, lwd = 2)
points(d$CO2_ppm[-(1:3)], post$BUGSoutput$median$ca, pch = 21, bg = "grey50",
       cex = 2, lwd = 2)

library(viridisLite)
cols = viridis(3)

## Plot showing sample and level distributions
dens = apply(post$BUGSoutput$sims.list$ca, 2, density)
dens.cl = apply(post.cl$BUGSoutput$sims.list$ca, 2, density)

ymax = 0
for(i in seq_along(dens.cl)){
  ymax = max(ymax, dens.cl[[i]]$y)
}
for(i in seq_along(dens)){
  ymax = max(ymax, dens[[i]]$y)
}

plot(0, 0, type = "n", xlim = c(100, 8000), ylim = c(0, ymax), 
     xlab = "ca", ylab = "Density")
for(i in seq_along(dens)){
  lines(dens[[i]], col = cols[data$level[i]])
}

for(i in seq_along(dens.cl)){
  lines(dens.cl[[i]], col = cols[i], lwd = 4)
}

legend("topright", legend = c(261, 292, 417), col = cols, lwd = 3, bty = "n",
       title = "Level")
