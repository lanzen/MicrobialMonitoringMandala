 ###################################################################################
 # Functions to handle ecological status predictions and ecogroups for a range of  #
 # indices, including custom abiotic BIs. Function to plot ML or QRS predictions   #
 # and calculate agreement (linear correlation and kappa)                          #
 #                                                                                 # 
 # Anders Lanzen Oct 2019                                                          #
 ###################################################################################

source("R/utils/ec_and_plot.R")

require(quantreg)
require(splines)
require(irr)

splinePredict = function(goodTITANInds, t, otus.train, md.train, imageOutDir = NA) {
  
  groups = data.frame(row.names=names(otus.train),value=rep(0,dim(otus.train)[2]))
  
  for (ti in c(1:dim(goodTITANInds)[1])){
    
    sp = row.names(goodTITANInds)[ti]
    titanGrp = goodTITANInds$maxgrp[ti]
    if (t=="Redox" | t=="Redox5yAvg" | t=="Redox10yAvg") titanGrp = 3-titanGrp
    
    ## Model distribution using 4th degree polynomial splines based on 95% quantiles 
    ## (Andersson et al 2008 and Nigel, except not sure how to check AIC and compare to 3 and 5 df)
    mt = data.frame(ab=otus.train[,sp],t=md.train[,t])
    
    bsp <- rq(ab ~ bs(t,df=6),data=mt,tau=.9)
    st = seq(min(md.train[,t]),max(md.train[,t]),length.out=1000)
    pred<-predict(bsp,data.frame(t=st))
    maxT = st[pred==max(pred)]
    ecoGroup = getECFromPeak(maxT,bi=t)
    
    ## If conflicting totally with TITAN, try with 5 DF
    
    if (is.na(ecoGroup) | (titanGrp==1 & ecoGroup>3) | (titanGrp==2 & ecoGroup<3)) {
      
      bsp <- rq(ab ~ bs(t,df=5),data=mt,tau=.9)
      st = seq(min(md.train[,t]),max(md.train[,t]),length.out=1000)
      pred<-predict(bsp,data.frame(t=st))
      maxT = st[pred==max(pred)]
      ecoGroup = getECFromPeak(maxT,bi=t)
      
      ## If still conflicting, do not assign to eco group
      if (is.na(ecoGroup) |(titanGrp==1 & ecoGroup>3) | (titanGrp==2 & ecoGroup<3)){
        print(paste("Warning: conflict for",sp,"TITAN =",titanGrp,"but predicted EC =",ecoGroup))
        ecoGroup = 0 
      }
    }
    
    groups[sp,"value"] = ecoGroup
    
    ## Plot distribution and check indicator status manually if needed 
    
    if (!is.na(imageOutDir)){
      
      change=NA
      if (goodTITANInds$maxgrp[ti]==2) {
        change <- "+"
      } else {
        if (goodTITANInds$maxgrp[ti]==1) change <- "-"
      }
      
      png(paste(imageOutDir,"/",sp,"_",t,".png",sep=""),width=400,height=400)
      
      plot(otus.train[,sp]~md.train[,t],
           main=paste(sp,change,"@",round(goodTITANInds$ienv.cp[ti],1)),
           sub=paste("Ecogroup",ecoGroup,"max @",round(maxT,1)),
           xlab=t,ylab="Abundance")
      abline(v=goodTITANInds$ienv.cp[ti], col="blue")
      abline(v=maxT, col="red", lty=2)
      lines(st,pred,col="orange")
      dev.off()
    }
  }
  return (groups)
}
