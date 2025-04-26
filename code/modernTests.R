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
prepMod(data)

system.time({
  post = jags.parallel(data, inits, parms, file.path(tempdir(), "model.txt"), n.chains = 3,
                       n.iter = 2e4, n.burnin = 1e4, n.thin = 1e3)
})

data = parseFranks(d)
prepMod(data)

system.time({
  post = jags.parallel(data, inits, parms, file.path(tempdir(), "model.txt"), n.chains = 3,
                       n.iter = 2e4, n.burnin = 1e4, n.thin = 1e3)
})



# Modern data, multi-site, multi-species ----
d = read.xlsx("data/stomata-franks_modernTests.xlsx", sheet = 1, startRow = 3)

## First w/o collapsing, single specimen interpretations
data = parseFranks(d, FALSE)
prepMod(data)

system.time({
  post = jags.parallel(data, inits, parms, file.path(tempdir(), "model.txt"), n.chains = 3,
                       n.iter = 2e4, n.burnin = 1e4, n.thin = 1e3)
})

## Now w/ collapsing
data = parseFranks(d)
system.time({
  post.cl = jags.parallel(data, inits, parms, "code/models/forwardFranksMulti.R", n.chains = 3,
                          n.iter = 2e6, n.burnin = 1e4, n.thin = 1e3)
})

