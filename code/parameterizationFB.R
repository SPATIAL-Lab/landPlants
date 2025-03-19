# Read supp data from Franks and Beerling 2009
fb = read.csv("data/FB09_supp.csv")

# Convert units
## 1 / mm^2 to 1 / m^2
fb$Density_m2 = fb$Density_mm2 * 1e6

## um^2 to m^2
fb$Size_m2 = fb$Size_microns / 1e12

# Calculate GC length
fb$GCL_m = sqrt(fb$Size_m2 / pi) * 2

# Density plots
## Stomatal density
plot(density(sqrt(fb$Density_m2), na.rm = TRUE))
range(fb$Density_m2, na.rm = TRUE)

## Stomatal GC length
plot(density(fb$GCL_m, na.rm = TRUE))
range(fb$GCL_m, na.rm = TRUE)

## Stomatal area
plot(density(fb$Size_m2 * fb$Density_m2, na.rm = TRUE))
range(fb$Size_m2 * fb$Density_m2, na.rm = TRUE)

# Parameterized distributions
## Stomatal area
plot(density(fb$Size_m2 * fb$Density_m2, na.rm = TRUE, from = 0, to = 1),
     main = "Stomatal Area")
lines(density(rbeta(1e6, 1.5, 20)), col = 2)
legend("topright", legend = c("Data", "Model"), col = 1:2, lwd = 1)

## Stomatal density, note 1e7 transformation for sampling efficiency
plot(density(fb$Density_m2 / 1e7, na.rm = TRUE, from = 0), main = "Stomatal Density")
alpha = 2
lines(density(rgamma(1e6, alpha, alpha / 25)), col = 2)
legend("topright", legend = c("Data", "Model"), col = 1:2, lwd = 1)

## Stomatal GC length, note 1e6 transformation for sampling efficiency
plot(density(fb$GCL_m * 1e6, na.rm = TRUE, from = 0), main = "GC Length")
alpha = 2
lines(density(rgamma(1e6, alpha, alpha / 25)), col = 2)
legend("topright", legend = c("Data", "Model"), col = 1:2, lwd = 1)
