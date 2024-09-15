model {
  
  # Likelihood ----
  D13C.obs ~ dnorm(D13C, 1 / D13C.sd ^ 2)
  v.obs ~ dnorm(v, 1 / v.sd ^ 2)
  a_st.obs ~ dnorm(a_st, 1 / a_st.sd ^ 2)
  l_leaf.obs ~ dnorm(l_leaf, 1 / l_leaf.sd ^ 2)
  AEO ~ dnorm(A / E, 1 / 0.05 ^ 2)
  ###need to maximize A/E
  
  # Equations ----
  ## boundary layer thickness ----
  d_bl = 4e-3 * sqrt(l_leaf / v_wind)
  
  ## leaf conductance is a function of leaf/stomatal geometry (Konrad eq 4)
  ### units m / s
  g_max = v * a_st * D_co2 / (d_st + v * a_st * (d_bl + d_as * tau_as^2 / n_as))
  
  ## transpiration is a function of g and env, units mol / m^2 / s
  ## use limiting transpiration rate to get limiting g
  g_lim = Elim / (1.6 * (w_sat - w_a))
  
  ## omega is stomatal regulation, either water limited or not
  omega_rel ~ dbeta(2, 2)
  omega = omega_rel * min(1, g_lim / g_max)
  
  ## finally g
  g = g_max * omega
  
  ## E
  E = 1.6 * g * (w_sat - w_a)

  ## Ci from Konrad eq 6, units mol / m^3
  C_i = 1 / (2 * g) * (g * (C_a - K) - (q - R_d) + 
                         sqrt((g * (C_a - K) - (q - R_d))^2 + 
                                4 * g * (g * K * C_a + q * Gam + K * R_d)))
  
  ## A from Ci and g, units mol / m^2 / s
  A = g * (C_a - C_i)
  
  ## D13C from Ci and Ca
  D13C = 4.4 + (27 - 4.4) * C_i / C_a
  
  # Intermediates ----
  p_sat = 6.1094 * exp(17.625 * t_air / (t_air + 243.04)) * 100  # saturation vapor pressure in hPa 
  w_sat = p_sat / (8.314 * (t_air + 273.15))  # saturation mol / m^3, ideal gas law
  w_a = w_sat * rh_air
  C_a = 101325 / (8.314 * (t_air + 273.15)) * CO2_a * 1e-6  # Atmospheric CO2 in mol / m^3, ideal gas law
  
  # Parameters ----
  v ~ dnorm(1e8, 1 / 1e7 ^ 2)T(0,)  # stomatal density in 1 / m^2
  d_st ~ dnorm(33.8e-6, 1 / 0.5e-6 ^ 2)T(0,)  #stomatal depth in m
  D_co2 = 1.55e-5  # diffusivity of CO2, T-dependent in m^2/sec
  a_st ~ dnorm(12e-12, 1 / 1e-12 ^ 2)T(0,) # stomatal pore area in m^2
  l_leaf ~ dnorm(8e-2, 1 / 1e-2 ^ 2)T(0,) # leaf length in m 
  v_wind = 1 # wind velocity in m / s
  d_as = 217.7e-6  # assimilating tissue thickness in m
  tau_as = 1.571  # assimilating tissue tortuosity dimensionless
  n_as = 0.35  # assimilating tissue porosity dimensionless
  t_air = 25  # air temperature in deg C
  rh_air ~ dunif(0.6, 0.9)  # air relative humidity
  K = 13925e-6  # MM constant in mol / m^3
  q ~ dnorm(7.3e-6, 1 / 0.5e-6 ^ 2)  # Carboxylation limit in mol / m^2 / s
  R_d = 0.16e-6  # Mitochondrial respiration rate in mol / m^2 / s
  Gam = 2408e-6  # Compensation pt in mol / m^3
  CO2_a ~ dnorm(350, 1 / 100 ^ 2)T(100,)  # Atmospheric CO2 ppm
}