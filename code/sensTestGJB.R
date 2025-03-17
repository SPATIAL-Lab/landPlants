# Test the forward model across range of CO2 and 2nd parameter ----
CO2_a = rep(seq(250, 1500, length = 25), 4)
Omega = c(rep(0, 25), rep(0.2, 25), rep(0.4, 25), 
      rep(0.6, 25))

source("code/models/forwardModel.R")
tst = fm(CO2_a = CO2_a, Omega = Omega) 

tst$CO2 = CO2_a
tst$Omega = Omega

png("out/senstest.png", 7, 5, "in", res = 300)
layout(matrix(1:4, nrow = 2))
par(mar = c(4, 4, 1, 1))
cats = unique(tst[, 7])
for(i in 1:4){
  dsub = tst[tst[, 7] == cats[1],]
  plot(dsub[, 6], dsub[, 1 + i], xlab = names(dsub)[6], 
       ylab = names(dsub)[1 + i], type = "l", ylim = range(tst[, 1 + i]))
  for(j in 2:length(cats)){
    dsub = tst[tst[, 7] == cats[j],]
    lines(dsub[, 6], dsub[, 1 + i], xlab = names(dsub)[6], 
          ylab = names(dsub)[1 + i], col = j)
  }
  if(i == 3){
    legend("bottomright", legend = paste(names(tst)[7], "=", cats), col = seq_along(cats), 
         lty = 1, bty = "n", bg = NULL)
  }
}
dev.off()

# Test initial JAGS implementation ----
library(R2jags)

## Baseline scenario
inits = NULL
parms = c("v", "D13C", "CO2_a", "q", "rh_air", "a_st", "d_st", 
          "g", "A", "E", "omega", "w_sat", "C_a", "C_i")

data = list("v.obs" = 100e6, "v.sd" = 1e6, "D13C.obs" = 18, "D13C.sd" = 0.3,
            "a_st.obs" = 12.3e-12, "a_st.sd" = 0.1e-12, "l_leaf.obs" = 84e-3,
            "l_leaf.sd" = 1e-3, "Elim" = 1e-4, "AEO" = 1)
p1 = jags(data, inits, parms, "code/JAGSmodel.R")

## Lower stomatal density
data = list("v.obs" = 50e6, "v.sd" = 1e6, "D13C.obs" = 22, "D13C.sd" = 0.3,
            "a_st.obs" = 12.3e-12, "a_st.sd" = 0.1e-12, "l_leaf.obs" = 84e-3,
            "l_leaf.sd" = 1e-3, "Elim" = 1e-4, "AEO" = 1)
p2 = jags(data, inits, parms, "code/JAGSmodel.R")

## Higher D13C
data = list("v.obs" = 200e6, "v.sd" = 1e6, "D13C.obs" = 18, "D13C.sd" = 0.3,
            "a_st.obs" = 12.3e-12, "a_st.sd" = 0.1e-12, "l_leaf.obs" = 84e-3,
            "l_leaf.sd" = 1e-3, "Elim" = 1e-4, "AEO" = 1)
p3 = jags(data, inits, parms, "code/JAGSmodel.R")

## Higher Elim
data = list("v.obs" = 100e6, "v.sd" = 1e6, "D13C.obs" = 16, "D13C.sd" = 0.3,
            "a_st.obs" = 12.3e-12, "a_st.sd" = 0.1e-12, "l_leaf.obs" = 84e-3,
            "l_leaf.sd" = 1e-3, "Elim" = 2e-4, "AEO" = 1)
p4 = jags(data, inits, parms, "code/JAGSmodel.R")

## Compare posterior means
p = rbind("Baseline" = p1$BUGSoutput$summary[,"mean"],
          "Low v High D13C" = p2$BUGSoutput$summary[,"mean"],
          "High D13C" = p3$BUGSoutput$summary[,"mean"],
          "Hig Elim" = p4$BUGSoutput$summary[,"mean"])
p

