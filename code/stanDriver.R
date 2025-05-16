library(openxlsx)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("code/helpers.R")

d = read.xlsx("data/stomata-franks_zhang_2024_P1.0.xlsx", sheet = 1, startRow = 3)

data = parseStan(d[10, ], FALSE)

post = stan("code/models/stupid.stan", chains = 1, data = data, open_progress = FALSE)

post
