model {
  
  for(i in 1:length(d13Cp[, 1])){
    # Likelihood
    d13Cp[i, 1] ~ dnorm(d13C_m[i], 1 / d13Cp[i, 2] ^ 2)
    Dab[i, 1] ~ dnorm(D[i], 1 / Dab[i, 2] ^ 2)
    GCLab[i, 1] ~ dnorm(Pl[i] / s1_m[i], 1 / GCLab[i, 2] ^ 2)
    GCWab[i, 1] ~ dnorm(l[i] / s2_m[species[i]], 1 / GCWab[i, 2] ^ 2)

    # Geometry
    Pl[i] = sqrt(SA[i] / D[i] / pi) * 2 
    
    # Franks model
    d13C_m[i] = d13Ca_m[level[i]] - D13C[i]
    D13C[i] = a + (b[i] - a) * ci[i] / ca[level[i]]
    ci[i] = ca[level[i]] - A[i] / gcop[i] - 1 / meso.scale[species[i]]
    
    # Based on data I've seen A should have noise added
    A[i] = (-q.b[i] - sqrt(q.b[i] ^ 2 - 4 * q.a[i] * q.c[i])) / (2 * q.a[i])
    
    q.a[i] = 1 / gcop[i] * (Ci0_m[species[i]] - gamma[i])
    q.b[i] = gamma[i] * (-2 * A0_m[species[i]] / gcop[i] + ca[level[i]] - 1 / meso.scale[species[i]] - 
                           2 * Ci0_m[species[i]] + 2 * gamma[i]) + 
      Ci0_m[species[i]] * (-A0_m[species[i]] / gcop[i] - ca[level[i]] + 1 / meso.scale[species[i]])
    q.c[i] = A0_m[species[i]] * (gamma[i] * (2 * ca[level[i]] - 2 / meso.scale[species[i]] - 
                                      Ci0_m[species[i]] - 2 * gamma[i]) +
                          Ci0_m[species[i]] * (ca[level[i]] - 1 / meso.scale[species[i]]))
    
    gcop[i] = (1 / gcop.g[i] + 1 / gb_m[species[i]]) ^ -1
    gcop.g[i] = gcmax[i] * gc.scale[species[i]]
    gcmax[i] = (d.v * D[i] * amax[i]) / (l[i] + ((pi / 2) * sqrt(amax[i] / pi))) / 1.6
    amax[i] = SA[i] * amax.scale[species[i]] / D[i] 
    
    # Priors
    ## Individual
    ### Stomatal area
    SA[i] ~ dbeta(1.5, 20)
    
    ### Stomatal density
    D[i] = d[i] * 1e7
    d[i] ~ dgamma(2, 2 / 25)

    l[i] ~ dunif(1e-6, 1e-4)

    s1_m[i] ~ dgamma(s1[i, 1] * s1.beta[i], s1.beta[i])
    s1.beta[i] = s1[i, 1] / s1[i, 2] ^ 2
    
  }
  
  ## Taxon
  for(i in 1:length(gb[, 1])){
    s2_m[i] ~ dgamma(s2[i, 1] * s2.beta[i], s2.beta[i])
    s2.beta[i] = s2[i, 1] / s2[i, 2] ^ 2
    
    amax.scale[i] ~ dbeta(s3[i, 1] * amax.v[i], (1 - s3[i, 1]) * amax.v[i]) # aka s3
    amax.v[i] = (s3[i, 1] * (1 - s3[i, 1])) / s3[i, 2] ^ 2 - 1
    
    gc.scale[i] ~ dbeta(s4[i, 1] * gc.v[i], (1 - s4[i, 1]) * gc.v[i]) # aka s4
    gc.v[i] = (s4[i, 1] * (1 - s4[i, 1])) / s4[i, 2] ^ 2 - 1
    
    meso.scale[i] ~ dbeta(s5[i, 1] * meso.v[i], (1 - s5[i, 1]) * meso.v[i]) # aka s5
    meso.v[i] = (s5[i, 1] * (1 - s5[i, 1])) / s5[i, 2] ^ 2 - 1
    
    Ci0_m[i] ~ dgamma(Ci0[i, 1] * Ci0.beta[i], Ci0.beta[i])
    Ci0.beta[i] = Ci0[i, 1] / Ci0[i, 2] ^ 2
    
    A0_m[i] ~ dgamma(A0[i, 1] * A0.beta[i], A0.beta[i])
    A0.beta[i] = A0[i, 1] / A0[i, 2] ^ 2
    
    gb_m[i] ~ dgamma(gb[i, 1] * gb.beta[i], gb.beta[i])
    gb.beta[i] = gb[i, 1] / gb[i, 2] ^ 2
  }

  ## Locality
  for(i in 1:length(d13Ca[, 1])){
    ca[i] = ca.s[i] * 1e3
    ca.s[i] ~ dunif(0.1, 8)
    d13Ca_m[i] ~ dnorm(d13Ca[i, 1], 1 / d13Ca[i, 2] ^ 2)
  }
  
  # Constants
  pi = 3.14159265
  a = 4.4
  d.v = 0.000940096
}
