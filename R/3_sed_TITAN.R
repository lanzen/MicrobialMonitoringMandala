## Anders Lanzén, 2019-09-22
## Identifies reliable indicator taxa using TITAN from metabarcoding data
## and models using quantile regression splines

## ****************************************************************************************
## Several operations in 1_sed_basic_analysis.R must have been run previous to this script!
## ****************************************************************************************

source("R/utils/ec_and_plot.R")
source("R/utils/splinePred.R")
require(TITAN2)
require(vegan)

# The maximum number of taxa, present in at least 3 samples, to analyse in TITAN (by top abundance)
TaxaToUse = 1000 

# The threshold for TITAN indicator reliability and "purity" model quantile regression splines 
# (can be increased to 0.99 to improve fidelity but at the cost of loosing indicators)
titan_threshold = .95

mds$StationSeasonYear = as.factor(paste(mds$Station,mds$Season, mds$Year,sep=""))
mds$col_plot = as.numeric(mds$Estuary)

# ----- TITAN analysis of OTUs and taxa ----

for (t in trainOn) {
  
  included_samples = (!is.na(mds[,t]) & mds$Environment=="Estuarine" & mds$Nacid=="DNA")
  extra = "_DNA"
  
  md.train = mds[included_samples,]
  otus.train = otus.s.ra.f[included_samples,]
  taxa.train = taxa.ass.ra.f[included_samples,]
  
  otus.train.pa = decostand(otus.train, method="pa")
  otus.train = otus.train[,colSums(otus.train.pa)>3]
  otus.train = otus.train[,order(colSums(otus.train), decreasing=T)]
  otus.train = otus.train[,c(1:TaxaToUse)]
  taxa.train.pa = decostand(taxa.train, method="pa")
  taxa.train = taxa.train[,colSums(taxa.train.pa)>3]
  
  otus.titan = titan(md.train[,t],otus.train,nBoot=100,numPerm=100)
  taxa.titan = titan(md.train[,t],taxa.train,nBoot=100,numPerm=100)
  write.csv(otus.titan$sppmax,paste("TITAN2_Splines/TITAN2_results/OTUs_TITAN_all_",t,".csv", sep=""))
  write.csv(taxa.titan$sppmax,paste("TITAN2_Splines/TITAN2_results/Taxa_TITAN_all_",t,".csv", sep=""))
}

# ----- Read TITAN results and analyse across training datasets ------

## Read output of TITAN, mark pure and reliable taxa

otuInd = matrix(ncol = length(trainOn), nrow = dim(otus.s.ra.f)[2], 
                dimnames=list(names(otus.s.ra.f),trainOn))

taxaInd = matrix(ncol = length(trainOn), nrow = dim(taxa.ass.ra.f)[2],
                 dimnames=list(names(taxa.ass.ra.f),trainOn))

otuGrp = matrix(ncol = length(trainOn), nrow = dim(otus.s.ra.f)[2], 
                dimnames=list(names(otus.s.ra.f),trainOn))

taxaGrp = matrix(ncol = length(trainOn), nrow = dim(taxa.ass.ra.f)[2],
                 dimnames=list(names(taxa.ass.ra.f),trainOn))

otuCP = matrix(ncol = length(trainOn), nrow = dim(otus.s.ra.f)[2], 
                dimnames=list(names(otus.s.ra.f),trainOn))

taxaCP = matrix(ncol = length(trainOn), nrow = dim(taxa.ass.ra.f)[2],
                 dimnames=list(names(taxa.ass.ra.f),trainOn))

for (t in trainOn) {
  tind = read.csv(paste("TITAN2_Splines/TITAN2_results/Taxa_TITAN_all_",t,".csv", sep=""),header=T,row.names=1)
  tindS = tind[tind$reliability>=titan_threshold & tind$purity>=titan_threshold & tind$obsiv.prob<.05,]
  for (ti in c(1:dim(tindS)[1])){
    taxaInd[row.names(tindS)[ti],t] = tindS$IndVal[ti]
    taxaGrp[row.names(tindS)[ti],t] = tindS$maxgrp[ti]
    taxaCP[row.names(tindS)[ti],t] = tindS$ienv.cp[ti]
  }
  
  oind = read.csv(paste("TITAN2_Splines/TITAN2_results/OTUs_TITAN_all_",t,".csv", sep=""),header=T,row.names=1)
  oindS = oind[oind$reliability>=titan_threshold & oind$purity>=titan_threshold & oind$obsiv.prob<.05,]
  for (ti in c(1:dim(oindS)[1])){
    otuInd[row.names(oindS)[ti],t] = oindS$IndVal[ti]
    otuGrp[row.names(oindS)[ti],t] = oindS$maxgrp[ti]
    otuCP[row.names(oindS)[ti],t] = oindS$ienv.cp[ti]
  }
}

## Check group membership for conflicts and check taxa presence in microgAMBI

mgAMBI = read.delim("~/projects/uGAMBI_RC/microgambi2017.tsv",sep="\t",row.names = 1)

taxaGrp = as.data.frame(taxaGrp)
taxaGrp$Consensus = 0
taxaGrp$Redox = 3-taxaGrp$Redox
taxaGrp$Redox5yAvg = 3-taxaGrp$Redox5yAvg
taxaGrp$Redox10yAvg = 3-taxaGrp$Redox10yAvg
taxaGrp$mgAMBI2017Group = NA
taxaGrp$mgAMBI2017Consensus = NA

for (rn in c(1:dim(taxaGrp)[1])){
  for (cn in c(1:dim(taxaGrp)[2])){
    grp = taxaGrp[rn,cn]
    if (!is.na(grp) & !is.na(taxaGrp$Consensus[rn])){
      if (taxaGrp$Consensus[rn]==0) {
        taxaGrp$Consensus[rn]=grp
      }
      else if (taxaGrp$Consensus[rn]!=grp) {
        taxaGrp$Consensus[rn] = NA
      }
    }
  }
  if ((row.names(taxaGrp)[rn]) %in% row.names(mgAMBI)){
    mGrp = mgAMBI[row.names(taxaGrp)[rn],]/6+1
    taxaGrp$mgAMBI2017Group[rn] = mGrp
    if (is.na(taxaGrp$Consensus[rn])) taxaGrp$mgAMBI2017Consensus[rn] = NA
    else taxaGrp$mgAMBI2017Consensus[rn] = (mGrp == taxaGrp$Consensus[rn])
  }
}
summary(taxaGrp$mgAMBI2017Consensus)
# Mode   FALSE    TRUE    NA's 
# logical      83     144     761 

otuGrp = as.data.frame(otuGrp)
otuGrp$Consensus = 0
otuGrp$Redox = 3-otuGrp$Redox
otuGrp$Redox5yAvg = 3-otuGrp$Redox5yAvg
otuGrp$Redox10yAvg = 3-otuGrp$Redox10yAvg
for (rn in c(1:dim(otuGrp)[1])){
  for (cn in c(1:dim(otuGrp)[2])){
    grp = otuGrp[rn,cn]
    if (!is.na(grp) & !is.na(otuGrp$Consensus[rn])){
      if (otuGrp$Consensus[rn]==0) {
        otuGrp$Consensus[rn]=grp
      }
      else if (otuGrp$Consensus[rn]!=grp) {
        otuGrp$Consensus[rn] = NA
      }
    }
  }
}


taxaInd[is.na(taxaInd)] = 0
taxaIndNorm = decostand(t(taxaInd), method="total")

otuInd[is.na(otuInd)] = 0
otuIndNorm = decostand(t(otuInd), method="total")

## ---- Spline prediction for taxa overview (and trivial, no cross-val performance eval) -----

## Read output of TITAN, mark pure and reliable taxa with p<0.05

# Groups for TITAN2 identified indicators, 
# reclassified as I -- V using quantile splines
otuGrp = matrix(ncol = length(trainOn), nrow = dim(otus.s.ra.f)[2], 
                dimnames=list(names(otus.s.ra.f),trainOn))

taxaGrp = matrix(ncol = length(trainOn), nrow = dim(taxa.ass.ra.f)[2], 
                 dimnames=list(names(taxa.ass.ra.f),trainOn))

for (t in trainOn) {
  print(paste("---",t,"---"))
  
  included_samples = (!is.na(mds[,t]) & mds$Environment=="Estuarine" & mds$Nacid=="DNA")
  extra = "_DNA"
  
  md.train = mds[included_samples,]
  otus.train = otus.s.ra.f[included_samples,]
  taxa.train = taxa.ass.ra.f[included_samples,]
  
  ## Read TITAN2 output and select pure and reliable indicators
  oind = read.csv(paste("TITAN2_Splines/TITAN2_results/OTUs_TITAN_all_",t,".csv", sep=""),header=T,row.names=1)
  oindS = oind[oind$reliability>=titan_threshold & oind$purity>=titan_threshold,] # & oind$obsiv.prob<.05,]
  
  tind = read.csv(paste("TITAN2_Splines/TITAN2_results/Taxa_TITAN_all_",t,".csv", sep=""),header=T,row.names=1)
  tindS = tind[tind$reliability>=titan_threshold & tind$purity>=titan_threshold,] # & oind$obsiv.prob<.05,]
  
  ## Iterate over all selected OTU indicators
  otuEGs = splinePredict(oindS,t, otus.train, md.train, 
                         imageOutDir="img/sed/Splines/OTUs/modelled")
  otuGrp[,t] = otuEGs$value
  taxaEGs = splinePredict(tindS,t, taxa.train, md.train, 
                          imageOutDir="img/sed/Splines/taxa/modelled")
  taxaGrp[,t] = taxaEGs$value
  
  ## Calculate BI values based on TITAN2 picked spline models
  
  otu_indicators = otuEGs[otuEGs$value>0,,drop=F] 
  taxa_indicators = taxaEGs[taxaEGs$value>0,,drop=F]
  
  otu_w_Ind = otus.train[,row.names(otu_indicators)]
  taxa_w_Ind = taxa.train[,row.names(taxa_indicators)]
  
  otu_BI_values = rep(0,dim(otus.train)[1])
  for (i in c(1:dim(otu_indicators)[1])){
    w = getWeight(otu_indicators[i,],bi=t)
    otu_BI_values = otu_BI_values + w*otu_w_Ind[,i]
  }
  otu_BI_values = otu_BI_values/rowSums(otu_w_Ind)
  
  taxa_BI_values = rep(0,dim(taxa.train)[1])
  for (i in c(1:dim(taxa_indicators)[1])){
    w = getWeight(taxa_indicators[i,],bi=t)
    taxa_BI_values = taxa_BI_values + w*taxa_w_Ind[,i]
  }
  taxa_BI_values = taxa_BI_values/rowSums(taxa_w_Ind)
  
  
  opreds = otu_BI_values[!is.na(otu_BI_values)]
  mdo.test = md.train[!is.na(otu_BI_values),]
  
  tpreds = taxa_BI_values[!is.na(taxa_BI_values)]
  mdt.test = md.train[!is.na(taxa_BI_values),]
  
  
  pdf(paste("img/sed/Splines/OTUs/",t,"_all_no_crossval_pred",extra,".pdf",sep=""),
      width=8,height=5)
  plot_ml(data = opreds, metadata = mdo.test, 
          xIndex = t, yIndex=t,
          aggreg = c("StationSeasonYear", "Estuary"), 
          title = "TITAN2+spline predicted BI values")
  dev.off()
  
  pdf(paste("img/sed/Splines/taxa/",t,"_all_no_crossval_pred",extra,".pdf",sep=""),
      width=8,height=5)
  plot_ml(data = tpreds, metadata = mdt.test, 
          xIndex = t, yIndex=t,
          aggreg = c("StationSeasonYear", "Estuary"), 
          title = "TITAN2+spline predicted BI values")
  dev.off()
}

# ----------  Write summarised results to files --------
# 
# tx <- otus.all[row.names(otuInd),"classification"]
# otuInd$taxonomy = tx

## Write TITAN2 summarised results
write.csv(otuInd, "TITAN2_Splines/OTUs_IndVal.csv")
write.csv(taxaInd, "TITAN2_Splines/Taxa_IndVal.csv")
write.csv(otuGrp, "TITAN2_Splines/OTUs_IndGroup.csv")
write.csv(taxaGrp, "TITAN2_Splines/Taxa_IndGroup.csv")
write.csv(otuCP, "TITAN2_Splines/OTUs_IndCP.csv")
write.csv(taxaCP, "TITAN2_Splines/Taxa_IndCP.csv")


# ------- Make heatmap for taxon importance = IndVal -------

require(gplots)

totalTaxaInd = colSums(taxaIndNorm)
taxaInd_chart = as.data.frame(taxaIndNorm[,order(totalTaxaInd,decreasing=TRUE)][,c(1:50)])
row.names(taxaInd_chart) = c("AMBI","PI","PI5y","Metals","HC","OM","Redox")
rowSums(taxaInd_chart)
# AMBI         PI       PI5y     Metals         HC         OM      Redox 
# 0.09853285 0.11076855 0.09960296 0.11932115 0.10890758 0.09682688 0.11537482 +

taxaInd_chart = decostand(taxaInd_chart,method="total")

# Get ecogroup as labels for heatmap, convert to roman numerals, and, importantly set IndVal to 0
# where ecogroup is zero, i.e. not pure or reliable or disagrees with TITAN
taxaInd_labels = taxaGrp[names(taxaInd_chart),]
for (i in c(1:50)){
  taxaInd_labels[i,] = as.character(as.roman(taxaInd_labels[i,]))
  taxaInd_chart[is.na(taxaInd_labels[i,]),i] = NA
}
min(taxaInd_chart,na.rm=T) #.0015
max(taxaInd_chart,na.rm=T) #.0024


pal <- colorRampPalette(c("black", "blue", "purple","red","yellow"))(n = 256)

pdf(paste("img/sed/TITAN2/Taxon_importance_heatmap",extra,".pdf",sep=""),
    width=6,height=10)
# heatmap(as.matrix(taxaInd_chart),col=pal,Colv = NA,
#         Rowv = NA,revC=TRUE,scale="col",cexRow=.8,cexCol=.8,
#         margins = c(4,10))
heatmap.2(as.matrix(sqrt(t(taxaInd_chart))),col=pal,dendrogram="none",Colv = NA,
          Rowv = NA,revC=F,scale="none",cexRow=.7,cexCol=.8,
          margins = c(5,10),
          density.info="none",trace="none",
          cellnote = taxaInd_labels,
          notecol="white",key=T,notecex=.6,
          keysize=2.2,
          key.par=list(mar=c(6,1,12,1)),
          na.color = "black"
)

dev.off()


# --------- Make heatmap for OTU importance = IndVal --------
totalotuInd = colSums(otuIndNorm)
otuInd = as.data.frame(otuInd)

otuInd_chart = as.data.frame(otuInd[order(totalotuInd,decreasing=TRUE),][c(1:50),])
names(otuInd_chart) = c("AMBI","PI","PI5y","Metals","HC","OM","Redox")
otuInd_labels = otuGrp[row.names(otuInd_chart),]
for (i in c(1:50)){
  otuInd_labels[i,] = as.character(as.roman(otuInd_labels[i,]))
  otuInd_chart[i,is.na(otuInd_labels[i,])] = NA
}

titanOTUs = row.names(otuInd_chart)

tx <- otus.all[row.names(otuInd_chart),"classification"]
bestTx = array(dim=length(tx))
for (i in 1:length(tx)) bestTx[i] <- tail(unlist(strsplit(as.character(tx[i]), split=";", fixed=TRUE)), 1)
row.names(otuInd_chart) = paste(gsub("SWARM_","",row.names(otuInd_chart)) ,bestTx)

pdf(paste("img/sed/TITAN2/OTU_importance_heatmap",extra,".pdf",sep=""),
    width=6.2,height=10)
heatmap.2(sqrt(as.matrix(otuInd_chart)),col=pal,dendrogram="none",Colv = NA,
          Rowv = NA,revC=F,scale="none",cexRow=.7,cexCol=.8,
          margins = c(5,10),
          density.info="none",trace="none",
          cellnote = otuInd_labels,
          notecol="white",key=T,notecex=.6,
          keysize=2.2,
          key.par=list(mar=c(6,1,12,1)),
          na.color="black"
)
dev.off()
