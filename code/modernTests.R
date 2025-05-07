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

# For samples with s3 = 1 & es3 = 0.05, sub in viable distribution parameters
mu = 0.97
sigma = 0.05
s3.v = (mu * (1 - mu)) / sigma ^ 2 - 1
s3.synth = rbeta(1e6, mu * s3.v, (1 - mu) * s3.v)
plot(density(s3.synth, from = 0, to = 1))
mean(s3.synth)
sd(s3.synth)
d$s3[d$s3 == 1] = mu

# Subset greenhouse and field data
d.green.500 = d[d$Stratigraphic.Level == 500,]
d.green.1k = d[d$Stratigraphic.Level == 1000,]
d.field = d[d$Stratigraphic.Level == 400,]

## Greenhouse single specimen interpretations
data.500 = parseFranks(d.green.500)
prepMod(data.500)

system.time({
  post.500 = jags.parallel(data, inits, parms, file.path(tempdir(), "model.txt"), 
                       n.chains = 3, n.iter = 2e4, n.burnin = 1e4)
})

s.500 = grep("500", d.green$Stratigraphic.Level)
s.1k = grep("1000", d.green$Stratigraphic.Level)
cas = c(500, 1000)
gen = unique(d.green$Genus)

pd = list()
for(i in seq_along(data$d13Cp[, 1])){
  pd[[i]] = density(post$BUGSoutput$sims.list$ca[, i])
}

pd.max = numeric(length(pd))
for(i in seq_along(data$d13Cp[, 1])){
  pd.max[i] = max(pd[[i]]$y)
}

plot(0, 0, type = "n", xlim = c(300, 3000), ylim = c(0, max(pd.max)),
     xlab = "ca (ppm)", ylab = "Density")
for(i in seq_along(data$d13Cp[, 1])){
  lines(pd[[i]], col = match(d.green$Genus[i], gen), 
        lwd = match(d.green$Stratigraphic.Level[i], cas))
}

plot(density(post$BUGSoutput$sims.list$ca[, 1]))



## Collapsing at level of CO2 treatment
data = parseFranks(d.green)
system.time({
  post = jags.parallel(data, inits, parms, file.path(tempdir(), "model.txt"), 
                       n.chains = 3, n.iter = 2e4, n.burnin = 1e4)
})

plot(0, 0, type = "n", xlim = c(300, 3000), ylim = c(0, 0.2),
     xlab = "ca (ppm)", ylab = "Density")
lines(density(post$BUGSoutput$sims.list$ca[, 1]))
lines(density(post$BUGSoutput$sims.list$ca[, 2]))
