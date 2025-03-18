parseFranks = function(d, condense = TRUE){
  ## Parse values from sheet, obs, parameters, constants

  ## Drop aggregate estimate rows based on missing stomatal density
  d = d[!is.na(d$Dab), ]
  
  data = list()
  
  data.names = c("d13Cp", "Dab", "GCLab", "GCWab")
  data.sd = c(0.2, 1.5e5, 3e-7, 1.5e-7)
  
  for(i in seq_along(data.names)){
    ci = match(data.names[i], names(d))
    d.sub = d[, ci:(ci + 1)]
    d.sub[is.na(d.sub[, 2]), 2] = data.sd[i]
    d.sub[d.sub[, 2] == 0, 2] = data.sd[i]
    data[[i]] = d.sub
  }
  
  mp.names = c("d13Ca", "A0", "CiCa0", "gb", "s1", "s2", "s3",
               "s4", "s5")
  mp.sd = c(0.5, 0.25, 0.05, 0.05, 0.05, 0.05, 0.025, 0.01, 0.001)
  
  for(i in seq_along(mp.names)){
    ci = match(mp.names[i], names(d))
    d.sub = d[, ci:(ci + 1)]
    d.sub[is.na(d.sub[, 2]), 2] = mp.sd[i]
    if(mp.names[i] == "s1"){
      # Special case if Pl is measured directly 
      d.sub[d.sub[, 2] == 0, 2] = 0.001
    } else{
      d.sub[d.sub[, 2] == 0, 2] = mp.sd[i]
    }
    data[[i + 4]] = d.sub
  }
  
  c.names = c("CO2_0", "b", "gamma")
  for(i in seq_along(c.names)){
    ci = match(c.names[i], names(d))
    d.sub = d[, ci]
    data[[i + 13]] = d.sub
  }
  
  names(data) = c(data.names, mp.names, c.names)
  
  data$Ci0 = data$CiCa0 * data$CO2_0
  names(data$Ci0) = c("Ci0", "eCi0")
  data = data[!(names(data) %in% c("CiCa0", "CO2_0"))]
  
  if(condense){
    ## Sites index
    ci = match("Stratigraphic.Level", names(d))
    strats = unique(d[, ci])
    data$level = match(d[, ci], strats)
    
    ## Condense site level parameters
    ### First occurrence of each strat level
    fo = match(strats, d[, ci])
    data$d13Ca = data$d13Ca[fo, ]
    
    ## Taxa index
    gs = paste(d$Genus, d$Species)
    species = unique(gs)
    data$species = match(gs, species)
    
    ## Condense species level parameters
    ### First occurrence of each species
    fo = match(species, gs)
    data$gb = data$gb[fo, ]
    data$A0 = data$A0[fo, ]
    data$Ci0 = data$Ci0[fo, ]
    data$s5 = data$s5[fo, ]
    data$s4 = data$s4[fo, ]
    data$s3 = data$s3[fo, ]
    data$s2 = data$s2[fo, ]
  } else{
    data$species = data$level = seq_len(nrow(d))
  }
  
  return(data)
}

inits = function() {
  list("ca.s" = runif(length(d13Ca[, 1]), 0.5, 2))  
}  
