model {
  
  # Likelihood
  d13C[1] ~ dnorm(d13C_m, 1 / d13C[2] ^ 2)
  Dobs[1] ~ dnorm(D, 1 / Dobs[2] ^ 2)
  GCL[1] ~ dnorm(Pl, 1 / GCL[2] ^ 2)
  GCW[1] ~ dnorm(l, 1 / GCW[2] ^ 2)
  
  # Franks model
  d13C_m = d13Ca - D13C
  D13C = a + (b - a) * ci / ca
  ci = ca - A / gcop - 1 / meso.scale
  
  A = (-q.b - sqrt(q.b ^ 2 - 4 * q.a * q.c)) / (2 * q.a)
  
  q.a = 1 / gcop * (Ci0 - gamma)
  q.b = gamma * (-2 * A0 / gcop + ca - 1 / meso.scale - 2 * Ci0 + 2 * gamma) +
    Ci0 * (-A0 / gcop - ca + 1 / meso.scale)
  q.c = A0 * (gamma * (2 * ca - 2 / meso.scale - Ci0 - 2 * gamma) +
                Ci0 * (ca - 1 / meso.scale))
  
  gcop = (1 / gcop.g + 1 / g.b) ^ -1
  gcop.g = gcmax * gc.scale
  gcmax = (d.v * D * amax) / (l + ((pi / 2) * sqrt(amax / pi))) / 1.6
  amax = (pi * (Pl / 2) ^ 2) * amax.scale
  
  # Priors
  Pl ~ dgamma(15, 1e6) # Pore length
  l ~ dgamma(15, 1e6) # Pore depth
  amax.scale ~ dbeta(100, 60)
  D ~ dgamma(7e1, 1e-6)
  gc.scale ~ dbeta(25, 100)
  ca ~ dunif(100, 2000)
  meso.scale ~ dbeta(130, 1e4)
  Ci0 ~ dgamma(250, 1)
  A0 ~ dgamma(5.8e3, 1e3)
  d13Ca ~ dnorm(-8.5, 1 / 0.1 ^ 2)
  g.b ~ dgamma(2e3, 1e3)
  
  # Constants
  d.v  = 0.000940096  # ratio of diffusivity of water vapor in air to the molar volume of air (mol m-1 s-1)
  gamma = 40  # CO2 compensation point (ppm) 
  a = 4.4  # discrimination against 13C during diffusion through air (per mil)
  b = 30  # discrimination against 13C due to carboxylation, mainly due to Rubisco (per mil) - make variable?
  pi = 3.14159265
}