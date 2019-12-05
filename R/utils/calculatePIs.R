###################################################################################
# Calculates abiotic pressure indices based on Aylagas et al (2017), but without  #
# including redox potential. Limits taken from Basque regional values (Menchaca   #
# 2014)                                                                           #
#                                                                                 # 
# Anders Lanzen 2019-11-12                                                        #
###################################################################################

calculatePIs = function(metadata){
  
  metadata$pi = 0
  metadata$piMetals = 0
  metadata$piHC = 0
  metadata$piOM = NA
  metalLimits = data.frame(row.names=c("Zn","Pb","Hg","Cd", "Cr", "Cu", "Ni"),
                           limits = c(249,78,.53,1,39,55,23))
  hcLimits = data.frame(row.names=c("TotalPAH","TotalPCB"), 
                        limits = c(1607,24.6))
  omLimit = 2.0
  allLimits = cbind(t(metalLimits), t(hcLimits), data.frame(OM=omLimit))
  
  geo = c(0,.5,1,2,4,8)
  for (s in row.names(metadata)){
    metalData = metadata[s,row.names(metalLimits)[row.names(metalLimits) %in% names(metadata)]]
    hcData = metadata[s,row.names(hcLimits)[row.names(hcLimits) %in% names(metadata)]]
    allPIData = cbind(metalData, hcData, data.frame(OM=metadata[s,"OM"]))
    
    mLim = metalLimits[row.names(metalLimits) %in% names(metadata),]
    hcLim = hcLimits[row.names(hcLimits) %in% names(metadata),]
    aLim = allLimits[,names(allLimits) %in% names(metadata)]
    
    hasPI = !is.na(allPIData)
    hasMData = !is.na(metalData)
    hasHCData = !is.na(hcData)
    
    # Get k-values (1/number of parameters)
    ka = 1/sum(hasPI)
    km = 1/sum(hasMData)
    khc = 1/sum(hasHCData)
    
    # Assign NA if no data
    if (sum(hasPI)==0) metadata[s,"pi"] = NA
    if (sum(hasMData)==0) metadata[s,"piMetals"] = NA
    if (sum(hasHCData)==0) metadata[s,"piHC"] = NA
    
    # Calculate total PI
    for (i in c(1:length(allPIData))[hasPI]){
      v = allPIData[,i]
      l = aLim[,i]
      for (j in c(1:5)){
        #print(paste(v,"<",l*geo[j+1],"?"))
        if (v<l*geo[j+1]){
          metadata[s,"pi"] = metadata[s,"pi"] + ka*(v-l*geo[j])/(l*geo[j+1] - l*geo[j])
          break
        } else {
          metadata[s,"pi"] = metadata[s,"pi"] + ka
        }
      }
    }
      
    # Calculate metal PI
    for (i in c(1:length(metalData))[hasMData]){
      v = metalData[,i]
      l = mLim[i]
      for (j in c(1:5)){
        if (v<l*geo[j+1]){
          metadata[s,"piMetals"] = metadata[s,"piMetals"] + km*(v-l*geo[j])/(l*geo[j+1] - l*geo[j])
          break
        } else {
          metadata[s,"piMetals"] = metadata[s,"piMetals"] + km
        }
      }
    }
      
    # Calculate HC PI
    for (i in c(1:length(hcData))[hasHCData]){
      v = hcData[,i]
      l = hcLim[i]
      for (j in c(1:5)){
        if (v<l*geo[j+1]){
          metadata[s,"piHC"] = metadata[s,"piHC"] + khc*(v-l*geo[j])/(l*geo[j+1] - l*geo[j])
          break
        } else {
          metadata[s,"piHC"] = metadata[s,"piHC"] + khc
        }
      }
    }
    
    # Calculate OM PI
    if (!is.na(metadata[s,"OM"])){
      om = metadata[s,"OM"]
      j=1
      piOM = 0
      while(om>omLimit*geo[j] & j<6){
        piOM = piOM + min(1,(om-omLimit*geo[j])/(omLimit*geo[j+1] - omLimit*geo[j]))  
        j=j+1
      }
      metadata[s,"piOM"] = piOM
    }
  }
  return(metadata)
}
