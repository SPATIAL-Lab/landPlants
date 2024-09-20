parseFranks = function(d){
  ## Parse values from sheet, obs, parameters, constants
  ## Need indexing on sample and on taxon 
  
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
  mp.sd = c(0.5, 0.25, 0.05, 0.05, 0.05, 0.025, 0.01, 0.001)
  
  for(i in seq_along(mp.names)){
    ci = match(mp.names[i], names(d))
    d.sub = d[, ci:(ci + 1)]
    d.sub[is.na(d.sub[, 2]), 2] = mp.sd[i]
    d.sub[d.sub[, 2] == 0, 2] = mp.sd[i]
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
  data = data[!(names(data) %in% c("CiCa0", "CO2_0"))]
  
  return(data)
}