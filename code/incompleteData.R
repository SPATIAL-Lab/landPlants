library(R2jags)
library(openxlsx)
source("code/helpers.R")
source("code/models/testModels.R")

# To save
parms = c("Pl", "l", "amax.scale", "D", "gc.scale", "ca", "meso.scale",
          "Ci0_m", "A0_m", "d13Ca_m", "A", "D13C", "gcop")

# Test dataset
d = read.xlsx("data/stomata-franks_zhang_2024_P1.0.xlsx", sheet = 1, startRow = 3)
data = parseFranks(d[4:8, ], FALSE)

# Run it
system.time({
  control = jags.parallel(data, NULL, parms, file.path(tempdir(), "fullFranks.txt"), 
                          n.chains = 4, n.iter = 2e6, n.burnin = 1e4, n.thin = 1e3)
})

# The reported uncertainty on d13Cp is unrealistic, try something more reasonable
data$d13Cp[, 2] = rep(0.3)

system.time({
  real.sigmad13C = jags.parallel(data, NULL, parms, file.path(tempdir(), "fullFranks.txt"), 
                          n.chains = 4, n.iter = 2e6, n.burnin = 1e4, n.thin = 1e3)
})

# Better convergence w/ the higher uncertainty
View(control$BUGSoutput$summary)
View(real.sigmad13C$BUGSoutput$summary)

# Very similar answers for all samples
for(i in seq_along(data$d13Cp[, 1])){
  plot(density(control$BUGSoutput$sims.list$ca[, i]), xlim = c(0, 8000),
       main = "", lwd = 2)
  lines(density(real.sigmad13C$BUGSoutput$sims.list$ca[, i]), col = 2, lwd = 2)
  
}

# Stomatal density only
## Simulate distribution for unknown d13Cp, use 1 sigma of 4 per mil 
data$d13Cp[, 1] = data$d13Ca[, 1] - 19
data$d13Cp[, 2] = rep(4)

## Drop GCW and GCL
data = data[names(data) != "GCLab"]
data = data[names(data) != "GCWab"]

# Run it
system.time({
  D.only = jags.parallel(data, NULL, parms, file.path(tempdir(), "DonlyFranks.txt"), 
                          n.chains = 4, n.iter = 2e6, n.burnin = 1e4, n.thin = 1e3)
})

# d13C only
data = parseFranks(d[4:8, ], FALSE)


plot(density(control$BUGSoutput$sims.list$ca[, 1]))
lines(density(D.only$BUGSoutput$sims.list$ca[, 1]))

for(i in seq_along(data$d13Cp[, i])){
  plot(density(control$BUGSoutput$sims.list$ca[, i]), xlim = c(0, 8000),
       main = "", lwd = 2)
  lines(density(D.only$BUGSoutput$sims.list$ca[, i]), col = 2, lwd = 2)
  
}
