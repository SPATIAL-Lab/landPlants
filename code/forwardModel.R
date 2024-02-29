
fm = function(# Parameters ----
              v = 85 * 1e3^2,  # stomatal density in 1/m^2
              d_st = 33.8 / 1e6,  #stomatal depth in m
              D_co2 = 1.55e-5,  # diffusivity of CO2, T-dependent in m^2/sec
              a_st = 12.3 / 1e6^2, # stomatal pore area in m^2
              d_bl = 0.66 / 1e3,  # boundary layer thickness in m
              d_as = 217.7 / 1e6,  # assimilating tissue thickness in m
              tau_as = 1.571,  # assimilating tissue tortuosity dimensionless
              n_as = 0.35,  # assimilating tissue porosity dimensionless
              t_air = 25,  # air temperature in deg C
              rh_air = 0.8,  # air relative humidity
              K = 13925 / 1e6,  # MM constant in mol / m^3
              q = 7.3 / 1e6,  # Carboxylation limit in mol / m^2 / s
              R_d = 0.16 / 1e6,  # Mitochondrial respiration rate in mol / m^2 / s
              Gam = 2408 / 1e6,  # Compensation pt in mol / m^3
              C_a = 2.68e25 / 6.023e23 * 350e-6  # Atmospheric CO2 in mol / m^3
              ){
  # Intermediates ----
  p_sat = exp(-(40700 / 8.3145) * (1 / (273 + t_air) - 1 / 373))  # saturation vapor pressure
  w_sat = 2.68e25 / 6.023e23 * p_sat  # saturation mol / m^3, assuming ideal gas and STP
  w_a = w_sat * rh_air
  
  # Equations ----
  ## stomatal conductance is a function of leaf/stomatal geometry (Konrad eq 4)
  g = v * d_st * D_co2 / (d_st + v * a_st * (d_bl + d_as * tau_as^2 / n_as))
  
  ## transpiration is a function of g and env
  E = 1.6 * g * (w_sat - w_a) 
  
  ## Ci from Konrad eq 6 
  C_i = 1 / (2 * g) * (g * (C_a - K) - (q - R_d) + 
                         sqrt((g * (C_a - K) - (q - R_d))^2 + 
                                4 * g * (g * K * C_a + q * Gam + K * R_d)))
  
  ## A from Ci and g
  A = g * (C_a - C_i)
  
  ## D13C from Ci and Ca
  D13C = 4.4 + (27 - 4.4) * C_i / C_a

  # Return ----
  return(data.frame(g, E, C_i, A, D13C))
}

fm(rh_air = c(0.6, 0.7, 0.8))
