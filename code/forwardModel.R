# Function ----
fm = function(## Parameters
              v = 85 * 1e3^2,  # stomatal density in 1 / m^2
              d_st = 33.8e-6,  #stomatal depth in m
              D_co2 = 1.55e-5,  # diffusivity of CO2, T-dependent in m^2/sec
              a_st = 12.3 / 1e6^2, # stomatal pore area in m^2
              d_bl = 0.66e-3,  # boundary layer thickness in m
              d_as = 217.7e-6,  # assimilating tissue thickness in m
              tau_as = 1.571,  # assimilating tissue tortuosity dimensionless
              n_as = 0.35,  # assimilating tissue porosity dimensionless
              t_air = 25,  # air temperature in deg C
              rh_air = 0.8,  # air relative humidity
              K = 13925e-6,  # MM constant in mol / m^3
              q = 7.3e-6,  # Carboxylation limit in mol / m^2 / s
              R_d = 0.16e-6,  # Mitochondrial respiration rate in mol / m^2 / s
              Gam = 2408e-6,  # Compensation pt in mol / m^3
              CO2_a = 350,  # Atmospheric CO2 ppm
              Omega = 1
              ){
  # Intermediates ----
  p_sat = 6.1094 * exp(17.625 * t_air / (t_air + 243.04)) * 100  # saturation vapor pressure in hPa 
  w_sat = p_sat / (8.314 * (t_air + 273.15))  # saturation mol / m^3, ideal gas law
  w_a = w_sat * rh_air
  C_a = 101325 / (8.314 * (t_air + 273.15)) * CO2_a * 1e-6  # Atmospheric CO2 in mol / m^3, ideal gas law
  
  # Equations ----
  ## stomatal conductance is a function of leaf/stomatal geometry (Konrad eq 4)
  ### units m / s
  g = v * a_st * (1 - Omega) * D_co2 / (d_st + v * a_st * (d_bl + d_as * tau_as^2 / n_as))
  
  ## transpiration is a function of g and env, units mol / m^2 / s
  E = 1.6 * g * (w_sat - w_a) 
  
  ## Ci from Konrad eq 6, units mol / m^3
  C_i = 1 / (2 * g) * (g * (C_a - K) - (q - R_d) + 
                         sqrt((g * (C_a - K) - (q - R_d))^2 + 
                                4 * g * (g * K * C_a + q * Gam + K * R_d)))
  
  ## A from Ci and g, units mol / m^2 / s
  A = g * (C_a - C_i)
  
  ## D13C from Ci and Ca
  D13C = 4.4 + (27 - 4.4) * C_i / C_a

  # Return ----
  return(data.frame(g, E, C_i, A, D13C))
}

# Konrad ----
## dw = water vapor deficit
## wi_sat = internal specific humidity
## wa = ambient specific humidity
## g = stomatal conductance 
## Ca = atmospheric CO2
## K = MM constant
## Gamma =
## q =
## Rd = 
## lambda = 
## a = 
## nu
## ast =
## dst = 
## Dco2 = Diffusivity of CO2 in m^2/s
## Dwv = Diffusivity of water vapor in m^2/s
## dbl =
## das = 
## tas = 
## nas = 
## A = assimilation rate
## E = specific transpiration
## t = air temperature in K
## Rgas = Ideal gas constant in kJ/mol
## p_atm = air pressure in J/m^3
## u = constant, mol/m^3
## v = constant, unitless
## Vm_w = molar volume of water in m^3/mol
## rho_w = water density in J*s^2/m^5
## grav = gravitational acceleration in m/s^2

## Water
wa = wrel * wi_sat
dw = wi_sat - wa

## Rubisco limited case
q = Vcmax_25 * exp(26.35 - 65.33 * (kJ / mol) / (Rgas * t))
Rd = Rd_25 * exp(18.72 - 46.39*(kJ / mol) / (Rgas * t));
W = q * kappa * Ca / (kappa * Ca * K * (1 + pO / Ko))

## Some enzymatic stuff
Jmax = Jmax_25 * exp(17.57 - 43.54 * (kJ/mol)/(Rgas * t))
Ph_PSII = 0.352 + 0.022 * (t - 273) - 3.4 * 10^(-4) * (t - 273) ^ 2
Q2 = Q * al_l * Ph_PSII * be_p
Th_PSII = -4.154 + 0.018 * t - 0.0003700000000 * (t - 273) ^ 2 # t in K
Jp = (Q2 + Jmax - sqrt((Q2 + Jmax) ^ 2 - 4 * Th_PSII * Q2 * Jmax)) / 
  2 / Th_PSII

## ?
Kc = exp(38.05 - 79.43 * (kJ/mol) / (Rgas * t)) * (mumol / mol) * 
  (p_atm / Rgas / t) * kc
Ko = exp(20.30 - 36.38 * (kJ/mol) / (Rgas * t)) * mmol/mol * 
  (p_atm / Rgas / t) * ko
pO = 210 * mmol/mol * (p_atm / Rgas / t)
Ga0 = exp(19.02 - 37.83 * (kJ/mol) / (Rgas * t)) * (mumol/mol) *
  (p_atm / Rgas / t)
ga_G = 1
Gamma = ga_G * Ga0
K_c = Kc * (1 + pO / Ko)

## Conversion of radiative power I_l (in J/m^2/s) to Quantum yield Q (in mol/m^2/s)
Q = I_l / 2e5  # equation unclear, is this multiplication or division?

## Physical quantities
Dco2 = (t / 273.15) ^ (1.8) * 1.33 * 10^(-5)
Dwv = (t / 273.15) ^ (1.8) * 2.13 * 10^(-5)
a = 1.6
Rgas = 8.3143e-3
p_atm = 101325.0 #sea level pressure
u = 2.035e10
v = 5306
wi_sat = u / t * exp(-v / t)
Vm_w = 18.015e-6
rho_w = 10e3
grav = 9.81

g = -1 / (Ca + K) ^ 2 * (sqrt(((K + Gamma) * q * 
                                 (Rd * K + Ca * Rd - Ca * q + q * Gamma) /
                                 (-Ca - K + lambda * a * dw) / 
                                 lambda / a / dw)) * 
                           (2 * a * dw * lambda - Ca - K) + 
                           Ca * Rd + Rd * K + 2 * q * Gamma - 
                           Ca * q + q * K)

nu = 1 / ast * (dst * g) / (Dco2 - (dbl + das * tas ^ 2 / nas) * g)

A = -(((K + Gamma) * q * (Rd * K + Ca * Rd - Ca * q + q * Gamma) / 
         (-Ca - K + lambda * a * dw) / 
         lambda / a / dw) ^ (1/2) * a * dw * lambda + 
        Rd * K + Ca * Rd - Ca * q + q * Gamma) / (Ca + K)

E = a * g * dw

# Franks ----

## Constants
d.v  = 0.000940096  # ratio of diffusivity of water vapor in air to the molar volume of air (mol m-1 s-1)
gamma = 40  # CO2 compensation point (ppm) 
a = 4.4  # discrimination against 13C during diffusion through air (per mil)
b = 30  # discrimination against 13C due to carboxylation, mainly due to Rubisco (per mil) - make variable?

Pl = GCL * s1  # Pl = stomatal pore length (m)
l = GCW * s2  # l = stomatal depth (m)
amax = (pi * (Pl / 2) ^ 2) * s3  # amax = the area of a fully-open stomatal pore
gcmax = (d.v * D * amax) / (l + ((pi / 2)*sqrt(amax / pi))) / 1.6  # gcmax = maximum stomatal conductance to CO2; Note: gcmax multiplied by 1.6 gives maximum conductance to water vapor, gwmax
gcop = gcmax * s4  # total operating stomatal conductance to CO2

A = A0 * (((ci.ca * ca - gamma) * (CiCa0 * CO2_0 + 2 * gamma)) / 
            ((ci.ca * ca + 2 * gamma) * (CiCa0 * CO2_0 - gamma)))

D13Cap = (d13Ca - d13Cp) / (1 + d13Cp / 1000)  # D13Cap = the stable carbon isotopic fractionation between leaf and atmosphere (per mil)
ci.ca = (D13Cap - a) / (b - a)


