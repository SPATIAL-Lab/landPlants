# Code for estimating CO2 with the Franks et al. model (2014, Geophysical Research Letters, 'New constraints on atmospheric CO2 concentration for the Phanerozoic', manuscript number 2014GL060457).

# Code is modified by Dana Royer from the r-scripts written by Dan Maxbauer for Metasequoia (Maxbauer et al., 2014 Geology 42: 1027-1030). This version (v2) modifies three elements of the original version posted as supplementary information to the Franks et al (2014) paper: 
# 1) if a resampling of 1 is used, the mean values of the inputs are used (not a resample); this change makes the code execute faster for runs where uncertainty bands are not wanted.
# 2) a more robust formulation for calculating photosynthetic rate (A) is used that takes into account the present-day ci/ca ratio of leaves (CiCa0); this piece was written by Jennifer Kowalczyk; resulting CO2 estimates are revised upward by at most ~10%.
# 3) we (DLR, JBK, and Joseph Milligan) added the option to run the code without calculating photosynthetic rate (A), using only the first main equation (i.e., Eq. 1 in Franks et al., 2014) to calculate CO2. This option can be used if the photosynthetic rate is independently known, for example with extant plants. With this option, the user should set the variable fixed_A to "yes" and set the variable A0 to the known photosynthetic rate; the variables CO2_0 and CiCa0 will not be used, and can be left blank.
# Updates 2) and 3) are further described in Kowalczyk et al. (in prep).

# Some modifications of the code may be needed depending on the species; e.g., Metasequoia requires the extra input parameter "inner rectangular length" (Maxbauer et al., 2014).

# The file "Franks_model_input.csv" contains all key input parameters and can be of any row length. Two examples are included: the first illustrates a species with stomata on only the abaxial surface (hypostomatous; this is the most common case) and the second illustrates a species with stomata on both leaf surfaces (amphistomatous).

# If the goal of the uncertainty analysis is to quantify confidence in the median CO2 estimate, standard errors of the mean (s.e.m.) are the most appropriate values to include in the input file.

# For r-novices, running the code should be simple:
# 1: set working directory (e.g., desktop) (File->Change dir...);
# 2: put input file into this directory;
# 3: paste all code in this file into console; the two output files will be created in your directory (a summary file and a file with all resampled CO2 estimates). Code with one sample and 10000 simulations takes several seconds to run on a PC with a 3.1GHz processor and 8GB RAM.

########################################################
### Description of columns in Franks_model_input.csv ###
########################################################
#sample: sample name

#Dab: stomatal density (m^-2) on abaxial surface (average over stomatal and non-stomatal areas).
#eDab: error in Dab.

#Dad: stomatal density (m^-2) on adaxial surface (average over stomatal and non-stomatal areas). Enter zero if stomata are absent.
#eDad: error in Dad. Enter zero if absent.

#GCLab: guard cell length (m) on abaxial surface. Note: If stomatal pore length (Pl) can be directly measured, enter Pl here (and its error) and use a 1:1 scaling between GCL and Pl (s1=1) and an associated error of 0 (es1=0).
#eGCLab: error in GCLab.

#GCLad: guard cell length (m) on adaxial surface. Enter zero if stomata are absent. Note: If stomatal pore length (Pl) can be directly measured, enter Pl here (and its error) and use a 1:1 scaling between GCL and Pl (s1=1) and an associated error of 0 (es1=0).
#eGCLad: error in GCLad. Enter zero if absent.

#GCWab: single guard cell width (m) on abaxial surface. Note: If GCW cannot be measured, it can be estimated from GCL with the scalings presented in Table S2 of Franks et al (2014); note that this table relates guard cell *pair* width to GCL.
#eGCWab: error in GCWab.

#GCWad: single guard cell width (m) on adaxial surface. Enter zero if stomata are absent. Note: If GCW cannot be measured, it can be estimated from GCL with the scalings presented in Table S2 of Franks et al (2014); note that this table relates guard cell *pair* width to GCL.
#eGCWad: error in GCWad. Enter zero if absent.

#d13Cp: ratio of 13C/12C isotopes in leaf material, relative to that in the PDB standard (per mil).
#ed13Cp: error in d13Cp.

#d13Ca: ratio of 13C/12C isotopes in (paleo-)atmosphere air, relative to that in the PDB standard (per mil). For Cenozoic material, the analysis of Tipple et al. (2010, Paleoceanography, doi:10.1029/2009PA001851) is helpful.
#ed13Ca: error in d13Ca. For Cenozoic material, the analysis of Tipple et al. (2010) is helpful.

#CO2_0: atmospheric CO2 concentration associated with A0 (ppm) (e.g., present-day value); this variable is assumed to have no error.  If fixed_A is set to 'yes', this variable is not used.

#A0: photosynthetic rate at CO2_0 (umol/m2/s). This can be measured on a fossil's nearest living relative. Alternatively, see Franks et al. (2014) for some generic scalings.  If fixed_A is set to yes, the value of A0 should be set to the photosynthetic rate A.
#eA0: error in A0.

#CiCa0: present-day Ci/Ca value.  This can be estimated from gas exchange or carbon isotope measurements on a living relative, or a typical value can be used (e.g. 0.65).  If fixed_A is set to 'yes', this variable is not used.
#eCiCa0: error in CiCa0.  If fixed_A is set to 'yes', this variable is not used.

#gb: boundary layer conductance to CO2 (mol/m2/s). Franks et al. (2014) suggests a generic value of 2 for typical conditions.
#egb: error in gb.

#s1: scaling from guard cell length (GCL) to stomatal pore length (Pl). See Table S2 in Franks et al. (2014) for some generic scalings.  s1 is equivalent to the term alpha in Table S2.
#es1: error in s1.

#s2: scaling from single guard cell width (GCW) to stomatal depth (l). In the typical case where guard cells have a circular cross-section, this scaling = 1.
#es2 error in s2.

#s3: scaling from the area of a circle with the diameter of pore length to a_max (maximum area of stomatal pore). See Table S2 in Franks et al. (2014) for some generic scalings.  s3 is equivalent to beta in Table S2.
#es3: error in s3

#s4: scaling from maximum conductance to CO2 (gcmax) to operational conductance to CO2 (gcop). This can be measured on a fossil's nearest living relative. Alternatively, Franks et al. (2014) suggests a generic scaling of 0.2.  s4 is equivalent to zeta in Equation 2 and Table S1 of Franks et al. (2014).
#es4: error in s4.

#s5: scaling from photosynthetic rate (A) to mesophyll conductance to CO2 (gm). Franks et al. (2014) suggests a generic scaling of 0.013.
#es5: error in s5.

#fixed_A: If set to 'yes', the code assumes a fixed photosynthetic rate equivalent to the input rate (A0) and solves only for CO2, using the first main equation (i.e. Eq. 1 in Franks et al., 2014).  This option can be used if A is known independently.  Otherwise, the code jointly solves for CO2 and photosynthetic rate (A) by iteritively solving the two main equations given in Franks et al. (2014) and Kowalczyk et al. (in prep.).
#########


#STEP 1: set constants and create matrices
resampleN <- 100  #number of resamples in Monte Carlo simulation; enter "1" if no uncertainty bands are wanted
#a low and high percentile are included in the summary output file when resampleN > 1
low_percentile <- 2.5; low_percentile_label <- paste0(low_percentile,"_percentile")
high_percentile <- 97.5; high_percentile_label <- paste0(high_percentile,"_percentile")

input <- read.csv("data/Franks_model_input.csv") #read input file
sampleN <- length(input[,"Dab"])  #determine the number of samples (=rows)

#create matrices for the two output files
if (resampleN==1) {
  CO2.summary <- matrix(nrow=sampleN, ncol=5, dimnames=list(input$sample,c("%converged","ci_ca","A","gctot","CO2")))
} else {
CO2.summary <- matrix(nrow=sampleN, ncol=7, dimnames=list(input$sample,c("%converged","ci_ca","A","gctot","CO2",low_percentile_label,high_percentile_label)))
}
CO2.resamples <- matrix(nrow=resampleN, ncol=sampleN, dimnames=list(1:resampleN,input$sample))

#create matrices (resampleN x 1) for some of the calculated variables
gcop.gb <- matrix(nrow=resampleN, ncol=1) #stomatal + boundary layer conductance to CO2 (for both leaf surfaces)
gctot <- matrix(nrow=resampleN, ncol=1) #stomatal + boundary layer + mesophyll conductance to cO2 (for both leaf surfaces)
ci.ca <- matrix(nrow=resampleN, ncol=1) #ratio of intercellular to atmospheric CO2 concentration

#constant parameters in model (see Franks et al., 2014); no errors are assumed for these parameters
d.v  <- 0.000940096  #ratio of diffusivity of water vapor in air to the molar volume of air (mol m-1 s-1)
gamma <- 40 #CO2 compensation point (ppm) 
a <- 4.4		#discrimination against 13C during diffusion through air (per mil)
b <- 30		#discrimination against 13C due to carboxylation, mainly due to Rubisco  (per mil)

for (j in 1:sampleN)  #loop that handles the number of input samples
{
  
  #STEP 2: fill input matrices with resamples
  #when resampleN = 1, the mean input values are used (no resamples)
  if (resampleN==1){
    Dab <- input[j,"Dab"]
    GCLab <- input[j,"GCLab"]
    GCWab <- input[j,"GCWab"]
    if (input[j,"Dad"]>0) {
      Dad <- input[j,"Dad"]
      GCLad <- input[j,"GCLad"]
      GCWad <- input[j,"GCWad"]
    }
    d13Cp <- input[j,"d13Cp"]
    d13Ca <- input[j,"d13Ca"]
    A0 <- input[j,"A0"]
    CiCa0 <- input[j,"CiCa0"]
    gb <- input[j,"gb"]
    
    s1 <- input[j,"s1"]
    s2 <- input[j,"s2"]
    s3 <- input[j,"s3"]
    s4 <- input[j,"s4"]
    s5 <- input[j,"s5"]
  } else {
    #when resampleN > 1, generate random normal distributions for every input parameter (with an error term) using the mean value and companion error
    Dab <- rnorm(resampleN, input[j,"Dab"], input[j,"eDab"])
    GCLab <- rnorm(resampleN, input[j,"GCLab"], input[j,"eGCLab"])
    GCWab <- rnorm(resampleN, input[j,"GCWab"], input[j,"eGCWab"])
    if (input[j,"Dad"]>0) {
      Dad <- rnorm(resampleN, input[j,"Dad"], input[j,"eDad"])
      GCLad <- rnorm(resampleN, input[j,"GCLad"], input[j,"eGCLad"])
      GCWad <- rnorm(resampleN, input[j,"GCWad"], input[j,"eGCWad"])
    }
    d13Cp <- rnorm(resampleN, input[j,"d13Cp"], input[j,"ed13Cp"])
    d13Ca <- rnorm(resampleN, input[j,"d13Ca"], input[j,"ed13Ca"])
    A0 <- rnorm(resampleN, input[j,"A0"], input[j,"eA0"])
    CiCa0 <- rnorm(resampleN, input[j,"CiCa0"], input[j,"eCiCa0"])
    gb <- rnorm(resampleN, input[j,"gb"], input[j,"egb"])
    
    s1 <- rnorm(resampleN, input[j,"s1"], input[j,"es1"])
    s2 <- rnorm(resampleN, input[j,"s2"], input[j,"es2"])
    s3 <- rnorm(resampleN, input[j,"s3"], input[j,"es3"])
    s4 <- rnorm(resampleN, input[j,"s4"], input[j,"es4"])
    s5 <- rnorm(resampleN, input[j,"s5"], input[j,"es5"])
  }
  
  CO2_0 <- input[j,"CO2_0"] #this is the one input parameter with no assumed error
  fixed_A <-input[j,"fixed_A"]
  
  
  #STEP 3: calculate resampled gcop.gb (total operating stomatal and boundary layer conductance to CO2; mol m-2 s-1) (see Franks et al., 2014):
  #    Pl = GCL*s1, where Pl = stomatal pore length (m)
  #    l = GCW*s2, where l = stomatal depth (m)
  #    amax = (pi*(Pl/2)^2)*s3, where amax = the area of a fully-open stomatal pore
  #    gcmax = (d.v*D*amax)/(l+((pi/2)*sqrt(amax/pi)))/1.6, where gcmax = maximum stomatal conductance to CO2; Note: gcmax multiplied by 1.6 gives maximum conductance to water vapor, gwmax
  #    gcop = gcmax*s4, equals the total operating stomatal conductance to CO2
  #calculation of resampled ci.ca (ratio of intercellular to atmospheric CO2)
  #    D13Cap = (d13Ca-d13Cp)/(1+d13Cp/1000), where D13Cap = the stable carbon isotopic fractionation between leaf and atmosphere (per mil)
  #    ci.ca = (D13Cap-a)/(b-a)
  for (i in 1:resampleN){
    gcop_abaxial <- (d.v*Dab[i]*(pi*(GCLab[i]*s1[i]/2)^2)*s3[i])/(GCWab[i]*s2[i]+((pi/2)*sqrt((pi*(GCLab[i]*s1[i]/2)^2)*s3[i]/pi)))/1.6*s4[i]
    if (input[j,"Dad"]>0) { #for amphistomatous species
      gcop_adaxial <- (d.v*Dad[i]*(pi*(GCLad[i]*s1[i]/2)^2)*s3[i])/(GCWad[i]*s2[i]+((pi/2)*sqrt((pi*(GCLad[i]*s1[i]/2)^2)*s3[i]/pi)))/1.6*s4[i]
      gcop.gb[i] <- ((1/gcop_abaxial)+(1/gb[i]))^-1+((1/gcop_adaxial)+(1/gb[i]))^-1
    } else {  #for hypostomatous species
      gcop.gb[i] <- ((1/gcop_abaxial)+(1/gb[i]))^-1
    }
    ci.ca[i] <- (((d13Ca[i]-d13Cp[i])/(1+d13Cp[i]/1000))-a)/(b-a)
  } #end of i loop
  
  
  #the two functions that relate photosynthetic rate (A) and atm CO2 (ca); they are solved iteratively (see Franks et al., 2014) unless fixed_A is set to 'yes', in which case only f1 is used. 
  f1 <- function(A){
    ca <- A/(((1/gcop.gb[i])+(1/(s5[i]*A)))^-1*(1-ci.ca[i])) 
  }  
  f2 <- function(ca){
    A <- A0[i]*(((ci.ca[i]*ca-gamma)*(CiCa0[i]*CO2_0+2*gamma))/((ci.ca[i]*ca+2*gamma)*(CiCa0[i]*CO2_0-gamma)))
  }
  
  #initial conditions for the iterative calculations
  ca0 <-matrix(data=f1(A0), nrow=resampleN, ncol=1) #calculate ca0 from A0 (where ca0 = atmospheric CO2 concentration associated with A0)
if (fixed_A == 'yes') {A1 <- A0} else A1 <- matrix(data=f2(ca0), nrow=resampleN, ncol=1)  #calculate A from ca0  
ca1 <- matrix(data=10000, nrow=resampleN, ncol=1)  #set ca1 (x=10000 ppm) artificially high so that the Monte Carlo loop doesn't fortuituously converge in the first step
  ca2 <- ca0
  converge <- matrix(data=0, nrow=resampleN, ncol=1)  #zero is the value that signals a non-converged run
  
  
  #STEP 4: run Monte Carlo simulation of the iterative calculations
  i <- 1  #reset i to 1 so that the iteration will start
  for(i in 1:resampleN)
  {
    iteration.count <- 0
    while((abs(ca2[i]-ca1[i]))>0.001 & iteration.count<=100)  #this while loop handles the iteration; it will stop after 100 iterations, signalling non-convergeance 
    {
      ca1[i] <- ca2[i] 
      if (fixed_A == 'yes') {A1[i]<-A0[i]} else A1[i] <- f2(ca1[i]) #once convergeance has been reached, this is the final A value
    ca2[i] <- f1(A1[i]) #once convergeance has been reached, this is the final CO2 value
      iteration.count <- iteration.count+1
    } #end of while loop
    if (iteration.count<100)  {converge[i] <- 1}  #this signals successful convergeance (=1)
    gctot[i] <- ((1/gcop.gb[i])+(1/(s5[i]*A1[i])))^-1
  } #end of resampleN loop [i]
  
  
  #STEP 5: populate the output files
  CO2.summary[j,"%converged"] <- sum(converge/resampleN*100)  #this tells you the percentage of the resampleN computations that converged to a solution (solving for CO2 and A simultaneously); sometimes when estimated CO2 is low, the equations do not converge.
  CO2.summary[j,"ci_ca"] <- median(ci.ca, na.rm=TRUE)
  CO2.summary[j,"A"] <- median(A1, na.rm=TRUE)
  CO2.summary[j,"gctot"] <- median(gctot, na.rm=TRUE)
  CO2.summary[j,"CO2"] <- median(ca2, na.rm=TRUE)
  if (resampleN > 1)  {
  CO2.summary[j,low_percentile_label] <- quantile(ca2, probs=low_percentile/100, na.rm=TRUE)
  CO2.summary[j,high_percentile_label] <- quantile(ca2, probs=high_percentile/100, na.rm=TRUE)
  }
  CO2.resamples[,j] <- ca2
} #end of sampleN loop [j]

write.csv(CO2.summary, "Franks_CO2_summary.csv")
write.csv(CO2.resamples, "Franks_CO2_resamples.csv")

## Testing quadratic solution for A
Ci0 = CiCa0 * CO2_0
ca = ca2

### Quadratic terms
q.a = 1 / gcop.gb * (Ci0 - gamma)
#q.b = -(ca * (Ci0 - gamma) + 2 * gamma * (Ci0 - gamma) + 
#  A0 / gcop.gb * (Ci0 + 2 * gamma))
#q.c = A0 * (ca - gamma) * (Ci0 + 2 * gamma)
q.b = gamma * (-2 * A0 / gcop.gb + ca - 1 / s5 - 2 * Ci0 + 2 * gamma) +
  Ci0 * (-A0 / gcop.gb - ca + 1 / s5)
q.c = A0 * (gamma * (2 * ca - 2 / s5 - Ci0 - 2 * gamma) +
              Ci0 * (ca - 1 / s5))

A = (-q.b - sqrt(q.b ^ 2 - 4 * q.a * q.c)) / (2 * q.a)


## Scratch

gcop.gb = (d.v * D * (pi * (Pl / 2) ^ 2) * amax.scale) / 
  (l +((pi / 2) * sqrt((pi * (Pl / 2) ^ 2) * amax.scale / pi))) / 
  1.6 * gc.scale

gctot[i] <- ((1 / gcop.gb) + (1 / (g.b * A1[i])))^-1
