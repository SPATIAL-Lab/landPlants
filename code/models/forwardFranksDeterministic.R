# Test to demonstrate that the Franks forward model math checks out 

library(openxlsx)
library(R2jags)
source("code/helpers.R")

# Using Xiaoqing's data
d = read.xlsx("data/stomata-franks_zhang_2024_P1.0.xlsx", sheet = 1, startRow = 3)

# Deterministic version of the forward model
frank = function(data, ca){
  
  # Constants
  pi = 3.14159265
  a = 4.4
  d.v = 0.000940096
  
  # Parameters
  d13Ca_m = data$d13Ca[1, 1]
  amax.scale = data$s3[1, 1]
  gc.scale = data$s4[1, 1]
  meso.scale = data$s5[1, 1]
  Ci0_m = data$Ci0[1, 1]  
  A0_m = data$A0[1, 1]
  gb_m = data$gb[1, 1]
  s2_m = data$s2[1, 1]
  s1_m = data$s1[1, 1]
  Pl = data$GCLab[1, 1] / s1_m
  l = data$GCWab[1, 1] / s2_m
  D = data$Dab[1, 1] * 1e7
  s2_m = data$s2[1, 1]
  s1_m = data$s1[1, 1]
  gamma = data$gamma
  b = data$b
  
  amax = (pi * (Pl / 2) ^ 2) * amax.scale
  gcmax = (d.v * D * amax) / (l + ((pi / 2) * sqrt(amax / pi))) / 1.6
  gcop.g = gcmax * gc.scale
  gcop = (1 / gcop.g + 1 / gb_m) ^ -1
  
  q.a = 1 / gcop * (Ci0_m - gamma)
  q.b = gamma * (-2 * A0_m / gcop + ca - 1 / meso.scale - 
                   2 * Ci0_m + 2 * gamma) + Ci0_m * (-A0_m / gcop - ca + 1 / meso.scale)
  q.c = A0_m * (gamma * (2 * ca - 2 / meso.scale - Ci0_m - 2 * gamma) +
                  Ci0_m * (ca - 1 / meso.scale))
  
  A = (-q.b - sqrt(q.b ^ 2 - 4 * q.a * q.c)) / (2 * q.a)
  
  ci = ca - A / gcop - 1 / meso.scale
  D13C = a + (b - a) * ci / ca
  d13C_m = d13Ca_m - D13C
  
  return(d13C_m)
}

# Function to run a test
ftest = function(d, ca, s){
  data = parseFranks(d[s, ], FALSE)
  if(nrow(data$Dab) == 0){
    stop("You picked a row with aggregate data")
  }
  d13C = frank(data, ca)
  
  ## Plot modeled plant d13C vs Ca
  plot(ca, d13C, type = "l", main = paste("Row", s))
  
  ## Measured plant d13C
  abline(h = data$d13Cp[1, 1], col = "red")
  
  ## Inferred 95% CI from inverse analysis
  lines(c(d$`CO2.low.(ppm)`[s], d$`CO2.high.(ppm)`[s]), 
        rep(data$d13Cp[1, 1], 2), col = "red", lwd = 4)
  
  ## Add JAGS result as ca density
  system.time({
    post = jags.parallel(data, NULL, c("ca"), "code/forwardFranksMulti.R", n.chains = 3,
                         n.iter = 5.01e6, n.burnin = 1e4, n.thin = 1e2)
  })
  dens = density(post$BUGSoutput$sims.list$ca)
  dens$y = min(d13C) + diff(range(d13C)) * dens$y / max(dens$y)
  lines(dens, col = "blue", lwd = 2)
  lines(quantile(post$BUGSoutput$sims.list$ca, c(0.025, 0.975)), rep(min(d13C), 2), 
        col = "blue", lwd = 3)
}

# Range of atmospheric CO2 values to be evaluated
ca = seq(100, 8000, by = 100)

ftest(d, ca, 4)
ftest(d, ca, 5)


data = parseFranks(d[r, ], FALSE)
d13C = frank(data, ca)
## Plot modeled plant d13C vs Ca
plot(ca, d13C, type = "l")
## Measured plant d13C
abline(h = data$d13Cp[1, 1], col = "red")
## Inferred 95% CI from inverse analysis
lines(c(d$`CO2.low.(ppm)`[r], d$`CO2.high.(ppm)`[r]), 
     rep(data$d13Cp[1, 1], 2), col = "red", lwd = 4)



# Row 5
r = 5
data = parseFranks(d[r, ], FALSE)
d13C = frank(data, ca)
lines(ca, d13C)
abline(h = data$d13Cp[1, 1], col = "red")
lines(c(d$`CO2.low.(ppm)`[r], d$`CO2.high.(ppm)`[r]), 
      rep(data$d13Cp[1, 1], 2), col = "red", lwd = 4)

