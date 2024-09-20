library(R2jags)

# Test data ----
## Invert PSM using Dana's test Ginko data

input = read.csv("data/Franks_model_input.csv")[1,]

data = list("d13C" = c(input$d13Cp, input$ed13Cp),
            "Dobs" = c(input$Dab, input$eDab),
            "GCL" = c(input$GCLab, input$eGCLab),
            "GCW" = c(input$GCWab, input$eGCWab))

parms = c("Pl", "l", "amax.scale", "D", "gc.scale", "ca", "meso.scale",
          "Ci0", "A0", "d13Ca", "A", "D13C", "gcop")

post = jags(data, NULL, parms, "code/forwardFranks.R", n.iter = 5e4)

View(post$BUGSoutput$summary)
sl = post$BUGSoutput$sims.list

plot(density(runif(1e6, 100, 2000)), ylim = c(0, 3e-3), main = "CO2")
lines(density(sl$ca), col = "red")
abline(v = mean(sl$ca), col = "red")

plot(density(rgamma(1e6, 15, 1e6)), ylim = c(0, 2e6), main = "Pl")
lines(density(sl$Pl), col = "red")

## Do the same conditioned only on d13C

data = list("d13C" = c(input$d13Cp, input$ed13Cp))

parms = c("Pl", "l", "amax.scale", "D", "gc.scale", "ca", "meso.scale",
          "Ci0", "A0", "d13Ca", "A", "D13C", "gcop")

post2 = jags(data, NULL, parms, "code/forwardFranks13C.R", n.iter = 5e5,
             n.burnin = 1e5, n.thin = 10)

View(post2$BUGSoutput$summary)
sl = post2$BUGSoutput$sims.list

plot(density(runif(1e6, 100, 2000)), ylim = c(0, 3e-3), main = "CO2")
lines(density(sl$ca), col = "red")
abline(v = mean(sl$ca), col = "red")

plot(density(rgamma(1e6, 15, 1e6)), main = "Pl")
lines(density(sl$Pl), col = "red")

# Template data ----

library(openxlsx)
source("code/helpers.R")

## Lei file works but has missing uncertainties

d = read.xlsx("data/stomata-franks_lei_2018.xlsx", sheet = 1, startRow = 3)

data = parseFranks(d)

parms = c("Pl", "l", "amax.scale", "D", "gc.scale", "ca", "meso.scale",
          "Ci0_m", "A0_m", "d13Ca_m", "A", "D13C", "gcop")

post = jags(data, NULL, parms, "code/forwardFranksMulti.R", n.iter = 5e4)

View(post$BUGSoutput$summary)
plot(d$CO2_ppm, post$BUGSoutput$median$ca)
abline(0, 1)

## None of these work w/o further tinkering

### Both have text in the d13Cp column
d = read.xlsx("data/stomata-franks_du_2018.xlsx", sheet = 1, startRow = 3)
d = read.xlsx("data/stomata-franks_li_2019.xlsx", sheet = 1, startRow = 3)

### Fixed A
d = read.xlsx("data/stomata-franks_montanez_2016.xlsx", sheet = 1, startRow = 3)

### Both missing CiCa0
d = read.xlsx("data/stomata-franks_zhou_2020.xlsx", sheet = 1, startRow = 3)
d = read.xlsx("data/stomata-franks_richey_2018.xlsx", sheet = 1, startRow = 3)
