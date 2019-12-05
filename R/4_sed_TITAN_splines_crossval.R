## Anders Lanzén, 2019-09-22
## Identifies reliable indicator taxa using TITAN from metabarcoding data
## and models using quantile regression splines, all using cross-validation
## by rotating estuaries

## ****************************************************************************************
## Several operations in 1_sed_basic_analysis.R must have been run previous to this script!
## ****************************************************************************************

source("R/utils/ec_and_plot.R")
source("R/utils/sml_compo.R")
source("R/utils/splinePred.R")

require(irr)
require(TITAN2)
require(quantreg)
require(splines)
require(vegan)

# The maximum number of taxa, present in at least 3 samples, to analyse in TITAN (by top abundance)
TaxaToUse = 1000 

# The threshold for TITAN indicator reliability and "purity" model quantile regression splines 
# (can be increased to 0.99 to improve fidelity but at the cost of loosing indicators)
titan_threshold = .95

mds$StationSeasonYear = as.factor(paste(mds$Station,mds$Season, mds$Year,sep=""))
mds$col_plot = as.numeric(mds$Estuary)

## Iterate over each predictor (AMBI, PI ...)

for (t in trainOn) {
  
  print(paste("*** Training on",t,"***"))
  
  included_samples = (!is.na(mds[,t]) & mds$Environment=="Estuarine" & mds$Nacid=="DNA")
  extra = "_DNA"
  
  md.trainAll = mds[included_samples,]
  # uncomment to treat azoic samples as AMBI 6 effectively limiting to range 0--6 instead of 0--7
  # md.train[md.train$AMBI==7,]$AMBI <- 6 
  otus.trainAll = otus.s.ra.f[included_samples,]
  taxa.trainAll = taxa.ass.ra.f[included_samples,]
  
  estuaries = unique(md.trainAll$Estuary)
  preds = data.frame(row.names=row.names(md.trainAll), 
                     byOTUs=rep(NA,dim(md.trainAll)[1]),
                     byTaxa=rep(NA,dim(md.trainAll)[1]))
  
  ## Leave-one-out cross-validation - OTUs:
  ## Iterate over each estuary, remove from training set and predict its value
  for (est in estuaries){
    
    print(est)
      
    ## Use the other estuaries to select indicator taxa w TITAN2
    md.train = md.trainAll[md.trainAll$Estuary!=est,]
    otus.train = otus.trainAll[md.trainAll$Estuary!=est,]
    otus.train.pa = decostand(otus.train, method="pa")
    otus.train = otus.train[,colSums(otus.train.pa)>5]
    otus.train = otus.train[,order(colSums(otus.train),decreasing=T)]
    otus.train = otus.train[,c(1:TaxaToUse)]
    
    otus.titan = titan(md.train[,t],otus.train,nBoot=100,numPerm=100) #ncpus=4,
    ot = as.data.frame(otus.titan$sppmax)
    
    ## Derive an AMBI-like index using regression splines with the picked indicators
    oindS = ot[ot$reliability>=titan_threshold & ot$purity>=titan_threshold,]
    otuEGs = splinePredict(oindS,t, otus.train, md.train)
    otu_indicators = otuEGs[otuEGs$value>0,,drop=F] 
    
    ## Calculate values for this speicific index for the left out estuary
    
    md.pred = md.trainAll[md.trainAll$Estuary==est,]
    otus.pred = otus.trainAll[md.trainAll$Estuary==est,]
    otu_w_Ind = otus.pred[,row.names(otu_indicators)]
    
    
    otu_BI_values = rep(0,dim(otus.pred)[1])
    for (i in c(1:dim(otu_indicators)[1])){
      w = getWeight(otu_indicators[i,],bi=t)
      otu_BI_values = otu_BI_values + w*otu_w_Ind[,i]
    }
    otu_BI_values = otu_BI_values/rowSums(otu_w_Ind)
    preds[md.trainAll$Estuary==est,"byOTUs"] = otu_BI_values
    
    # Write predictions to disk
    write.table(preds,paste("TITAN2_Splines/SplinePreds/",t,".tsv",sep=""), 
                sep="\t",quote=F, col.names=NA)
  }
  
  ## Compare the results based on impact indeces (leave-one-out validation)
  pdf(paste("img/sed/Splines/OTUs/",t,"_pred",extra,".pdf",sep=""),width=8,height=5)
  plot_ml(data = preds$byOTUs[!is.na(preds$byOTUs)], 
          metadata = md.trainAll[!is.na(preds$byOTUs),], 
          xIndex = t, yIndex=t,
          aggreg = c("StationSeasonYear", "Estuary"), 
          title = "Spline predicted BI values")
  dev.off()
  
  ## Leave-one-out cross-validation - taxa:
  for (est in estuaries){
    
    print(est)
    
    ## Use the other estuaries to select indicator taxa w TITAN2
    md.train = md.trainAll[md.trainAll$Estuary!=est,]
    taxa.train = taxa.trainAll[md.trainAll$Estuary!=est,]
    taxa.train.pa = decostand(taxa.train, method="pa")
    taxa.train = taxa.train[,colSums(taxa.train.pa)>10]
    
    taxa.titan = titan(md.train[,t],taxa.train,nBoot=100,numPerm=100,memory=T) #ncpus=4
    tt = as.data.frame(taxa.titan$sppmax)
    
    ## Derive an AMBI-like index using regression splines with the picked indicators
    tindS = tt[tt$reliability>=titan_threshold & tt$purity>=titan_threshold,]
    taxaEGs = splinePredict(tindS,t, taxa.train, md.train)
    taxa_indicators = taxaEGs[taxaEGs$value>0,,drop=F]
    
    ## Calculate values for this speicific index for the left out estuary
    
    md.pred = md.trainAll[md.trainAll$Estuary==est,]
    taxa.pred = taxa.trainAll[md.trainAll$Estuary==est,]
    taxa_w_Ind = taxa.pred[,row.names(taxa_indicators)]
    
    taxa_BI_values = rep(0,dim(taxa.pred)[1])
    for (i in c(1:dim(taxa_indicators)[1])){
      w = getWeight(taxa_indicators[i,],bi=t)
      taxa_BI_values = taxa_BI_values + w*taxa_w_Ind[,i]
    }
    taxa_BI_values = taxa_BI_values/rowSums(taxa_w_Ind)
    preds[md.trainAll$Estuary==est,"byTaxa"] = taxa_BI_values
    
  }
  
  
  pdf(paste("img/sed/Splines/taxa/",t,"_pred",extra,".pdf",sep=""),width=8,height=5)
  plot_ml(data = preds$byTaxa[!is.na(preds$byTaxa)], 
          metadata = md.trainAll[!is.na(preds$byTaxa),], 
          xIndex = t, yIndex=t,
          aggreg = c("StationSeasonYear", "Estuary"), 
          title = "Spline predicted BI values")
  dev.off()
  
  # Write predictions to disk
  write.table(preds,paste("TITAN2_Splines/SplinePreds/",t,".tsv",sep=""), 
              sep="\t",quote=F, col.names=NA)
  
}
  

