## Anders Lanzén, 2019-09-22
## Reads morphosp. data and uses same ML, TITAN2 and splines methods as other scripts
## Requires separate handling due to the way samples have to be picked

## ****************************************************************************************
## Several operations in 1_sed_basic_analysis.R must have been run previous to this script!
## ****************************************************************************************

source("R/utils/sml_compo.R")
source("R/utils/ec_and_plot.R")
source("R/utils/splinePred.R")
require(irr)
require(TITAN2)
require(quantreg)
require(splines)
require(irr)
require(vegan)

# Calculate mean value for all replcates from same station, season and year
mds.nr.dna$StationSeasonYear = as.factor(paste(mds.nr.dna$Station,
                                               mds.nr.dna$Season, 
                                               mds.nr.dna$Year,sep=""))
mds.nr.dna$col_plot = as.numeric(mds.nr.dna$Estuary)

## Read morpho data
morpho = read.csv("morphology/morpho_table_all.csv",header=T,row.names=1)
morpho = morpho[rowSums(morpho)>0,]

# ------- Mantel test 16S v morpho ------

## Limit morpho data to samples where metabarcoding data exist
mdm = mds.nr.dna[(row.names(mds.nr.dna) %in% row.names(morpho)),]
morpho.m = morpho[row.names(mdm),]
otus.m = otus.s.ra.f.nr.dna[row.names(mdm),]
table(row.names(otus.m) == row.names(morpho.m))

mantel(vegdist(otus.m),vegdist(morpho.m))
# Mantel statistic r: 0.3487 ***

mantel(vegdist(otus.m),vegdist(decostand(morpho.m,method="total")))
# Mantel statistic r: 0.3196 ***

mdm2 = mdm[!is.na(mdm$piMetals),]
morpho.m2 = morpho[row.names(mdm2),]
otus.m2 = otus.s.ra.f.nr.dna[row.names(mdm2),]
mantel(vegdist(otus.m2),vegdist(morpho.m2))

par.dist = dist(decostand(mdm2[,c("Mud","piOM","piMetals","piHC","WaterO2",
                                 "Redox5yAvg","PointDepthBM","SalClass")],method="range"))

mantel(vegdist(otus.m2),par.dist)
#Mantel statistic r: 0.4439 ***

mantel(vegdist(morpho.m2),par.dist)
#Mantel statistic r: 0.3626 ***
# --> Metabarcoding community dissimilarities corresponds better to 
#physicochemical parameters than morpho do

mantel.partial(vegdist(otus.m2),vegdist(morpho.m2),par.dist)
#Mantel statistic r: 0.2176  ***
# -> communities depend on each other or on variables not modelled
# Both morpho and prok. are more structured by parameters than each other, so in different dir.
# But values are too close to call.

# ------- Morpho NMDS ---------

nmds.ms = metaMDS(morpho.m,distance="bray",engine="isoMDS")

efm = envfit(nmds.ms, mdm[,c(16, 21,25:30,32:34,39:44,50,54,56,73:75)],na.rm=T)

mcolour = rb[as.numeric(mdm$Estuary)] 
mcolU = rb[sort(as.numeric(unique(mdm$Estuary)))]
pch=as.numeric(mdm$Environment)-1
pch[mdm$Season=="Su" & mdm$Nacid=="DNA" & mdm$Environment=="Estuarine"] = 5
pch[mdm$Season=="Su" & mdm$Nacid=="RNA"] = 18
pch[mdm$Season=="Wi" & mdm$Nacid=="RNA" & mdm$Environment=="Estuarine"] = 16
pch[mdm$Season=="Wi" & mdm$Nacid=="RNA" & mdm$Environment=="Coastal"] = 15

svg("img/sed/morpho_NMDS_w_envfit.svg")

ordiplot(nmds.ms,type="none")
plot(efm,p.max=.001,cex=.5,col="purple")
points(nmds.ms,display="sites",col=mcolour, pch=pch,cex=.7,lwd=1)
legend("bottomright",box.lwd=0.5, pch=c(1,5,0,2,3),
       legend=c("Estuarine Wi","Estuarine Su","Coastal","Offshore","Port"),
       ncol=1,cex=.6)
legend("topright",box.lwd=0.5, pch=c(1),legend=sort(unique(mdm$Estuary)),col=colU,
       ncol=1,cex=.6)

ordihull(nmds.ms,groups=mdm$Environment,show.groups=c("Coastal","Offshore"),
         col="darkgrey")


dev.off()


# ------ 1) TITAN2 and splines cross-validation ------
for (t in trainOn) {
  
  print(paste("*** Training on",t,"***"))
  
  included_samples = (!is.na(mds.nr.dna[,t]) & mds.nr.dna$Environment=="Estuarine")
  md.trainAll = mds.nr.dna[included_samples,]
  
  ## Check that morpho data is available for selected training
  row.names(md.trainAll) = gsub("2017A","2017",row.names(md.trainAll))
  row.names(md.trainAll) = gsub("_PS","",row.names(md.trainAll))
  table(row.names(md.trainAll) %in% row.names(morpho))
  row.names(md.trainAll)[!(row.names(md.trainAll) %in% row.names(morpho))]
  
  ## Make subset of morpho for ML training
  md.trainAll = md.trainAll[(row.names(md.trainAll) %in% row.names(morpho)),]
  morpho.trainAll = morpho[row.names(md.trainAll),]
  
  estuaries = unique(md.trainAll$Estuary)
  preds = data.frame(row.names=row.names(md.trainAll), 
                     bymorpho=rep(NA,dim(md.trainAll)[1]))
  
  ## Leave-one-out cross-validation:
  ## Iterate over each estuary, remove from training set and predict its value
  for (est in estuaries){
    
    print(est)
    
    ## Use the other estuaries to select indicator taxa w TITAN2
    md.train = md.trainAll[md.trainAll$Estuary!=est,]
    morpho.train = morpho.trainAll[md.trainAll$Estuary!=est,]
    
    morpho.train.pa = decostand(morpho.train, method="pa")
    morpho.train = morpho.train[,colSums(morpho.train.pa)>3]
    morpho.train = morpho.train[,order(colSums(morpho.train))]
    
    morpho.titan = titan(md.train[,t],morpho.train,nBoot=100,numPerm=100)
    mt = as.data.frame(morpho.titan$sppmax)
    
    ## Derive an AMBI-like index using regression splines with the picked indicators
    mIndS = mt[mt$reliability>=.95 & mt$purity>=.95,]
    
    morphoEGs = splinePredict(mIndS,t, morpho.train, md.train)
    morpho_indicators = morphoEGs[morphoEGs$value>0,,drop=F] 
    
    ## Calculate values for this speicific index for the left out estuary
    
    md.pred = md.trainAll[md.trainAll$Estuary==est,]
    morpho.pred = morpho.trainAll[md.trainAll$Estuary==est,]
    morpho_w_ind = morpho.pred[,row.names(morpho_indicators)]
    morpho_BI_values = rep(0,dim(morpho.pred)[1])
    
    for (i in c(1:dim(morpho_indicators)[1])){
      w = getWeight(morpho_indicators[i,],bi=t)
      morpho_BI_values = morpho_BI_values + w*morpho_w_ind[,i]
    }
    morpho_BI_values = morpho_BI_values/rowSums(morpho_w_ind)
    preds[md.trainAll$Estuary==est,"bymorpho"] = morpho_BI_values
    
  }
  
  ## Compare the results based on AMBI-like indeces (leave-one-out validation)
  pdf(paste("img/sed/Splines/morpho/",t,"_pred",extra,".pdf",sep=""),width=8,height=5)
  plot_ml(data = preds$bymorpho[!is.na(preds$bymorpho)], 
          metadata = md.trainAll[!is.na(preds$bymorpho),], 
          xIndex = t, yIndex=t,
          aggreg = c("StationSeasonYear", "Estuary"), 
          title = "Spline predicted BI values")
  dev.off()
  
  # Print predictions to disk
  write.table(preds,paste("TITAN2_Splines/SplinePreds/",t,"_morpho.tsv",sep=""), 
              sep="\t",quote=F, col.names=NA)
}



## ---- 2 MACHINE LEARNING (INCL. MEAN OF SPLINE and ML) -------

# morpho importance across training sets matrices 

morphoImpRF = matrix(ncol = length(trainOn), nrow = dim(morpho)[2],
                 dimnames=list(names(morpho),trainOn))

# Train ML classifier with cross-val and apply TITAN2 without, 
# across training parameters

for (t in trainOn){
  print("")
  print(paste("***",t,"***"))
  print("Running ML")
  
  included_samples = (!is.na(mds.nr.dna[,t]) & mds.nr.dna$Environment=="Estuarine")

  md.train = mds.nr.dna[included_samples,]
  
  ## Check that morpho data is available for selected training
  row.names(md.train) = gsub("2017A","2017",row.names(md.train))
  row.names(md.train) = gsub("_PS","",row.names(md.train))
  table(row.names(md.train) %in% row.names(morpho))
  row.names(md.train)[!(row.names(md.train) %in% row.names(morpho))]
  
  ## Make subset of morpho for ML training
  md.train = md.train[(row.names(md.train) %in% row.names(morpho)),]
  morpho.train = morpho[row.names(md.train),]
  dim(morpho.train)
  
  # ---- Predict RF classifier on all datasets with leave-one-out testing (on estuary name) ----
  
  preds <- sml_compo(otu_table = morpho.train, metadata = md.train, index=t,
                     cross_val = "Estuary", 
                     algo = "RF", optim_overfit = T)
  
  pdf(paste("img/sed/ML/morpho/",t,"_RF_pred",extra,".pdf",sep=""),width=8,height=5)
  status <- plot_ml(data = preds, metadata = md.train, xIndex = t, yIndex=t,
                    aggreg = c("StationSeasonYear", "Estuary"), 
                    title = "Morphotaxonomy (estuarine)")
  
  dev.off()
  
  ## Get Spline pred data and make average:
  
  sPred = read.table(paste("TITAN2_Splines/SplinePreds/",t,"_morpho.tsv",sep=""), 
              sep="\t",row.names=1,header=T)
  
  mdx2 = md.train[!is.na(sPred$bymorpho),]
  preds_w_spline = preds[!is.na(sPred$bymorpho)]
  sPred = sPred[!is.na(sPred$bymorpho),,drop=F]
  predx2 = rowMeans(data.frame(splines=sPred,rf=preds_w_spline))
  
  pdf(paste("img/sed/ML/morpho/",t,"_RF_Spline_special_pred",extra,".pdf",sep=""),
      width=8,height=5)
  
  status <- plot_ml(data = predx2, metadata = mdx2, xIndex = t, yIndex=t,
                    aggreg = c("StationSeasonYear", "Estuary"), 
                    title = "Morphotaxonomy (estuarine)")
  
  dev.off()
  
  ## RF Variable importance
  
  tt = md.train[,t]
  
  mod <- ranger(tt ~ ., data=morpho.train, 
                mtry=floor(dim(morpho.train)[2]/3), num.trees = 300, 
                importance= "impurity", write.forest = T)
  morphoImpRF[,t] = mod$variable.importance
  
  # Most 30 important for graphs
  imp <- tail(sort(mod$variable.importance), 30)
  
  
  pdf(paste("img/sed/ML/morpho/",t,"_MorphoTaxon_importance",extra,".pdf",sep=""),width=5,height=6)
  p <- barplot(imp, horiz = T, xlab="Variable importance", 
               main=paste("Taxon importance "), las=2, cex.names = 0.6)
  
  text(0, p, labels=names(imp), pos=4, cex=0.6)
  
  dev.off()
  
  # ---------- TITAN2 analysis (to get importance - not cross-val) ------
  
 
  morpho.train.pa = decostand(morpho.train, method="pa")
  morpho.train = morpho.train[,colSums(morpho.train.pa)>3]
  
  print ("Running TITAN...")
  morpho.titan = titan(md.train[,t],morpho.train,nBoot=100,numPerm=100)
  write.csv(morpho.titan$sppmax,paste("TITAN2_Splines/TITAN2_results/morpho_TITAN_all_",t,".csv", sep=""))
  
}
  

# --- 3) Look at important morphotaxa importance in RFs across training sets ----

## Important taxa in RF

# Set NA values to zero
morphoImpRF[is.na(morphoImpRF)] = 0
morphoImpRFNorm = decostand(as.data.frame(t(morphoImpRF)),method="total")
totalmorphoImpRF = colSums(morphoImpRFNorm)
morphoImpRF_chart = morphoImpRFNorm[,order(totalmorphoImpRF,decreasing=TRUE)][,c(1:50)]

rowSums(morphoImpRF_chart) 
# 80% of PI Metals to 92 % of Redox

#using a red and blue colour scheme without traces and scaling each row
library("RColorBrewer")
library("gplots")
pal <- colorRampPalette(c("black", "blue", "purple","red"))(n = 299)

pdf(paste("img/sed/ML/Morphotaxa_importance_heatmap",extra,".pdf",sep=""),
    width=5.5,height=9)
heatmap((t(as.matrix(morphoImpRF_chart))),col=pal,Colv = NA,
        Rowv = NA,revC=TRUE,scale="col",cexRow=.8,cexCol=.8,
        margins = c(0,10))

dev.off()


# --- Indicators by TITAN across training sets, predict a BI ------

## Read AMBI for reference
AMBI = read.delim("~/projects/AMBI/ambi_oct2018.tsv",sep="\t",header=F,col.names=c("taxon","group"))
AMBI = AMBI[!duplicated(AMBI$taxon),]
row.names(AMBI) = AMBI$taxon
AMBI$group[AMBI$group==0] = 1
AMBI$group[AMBI$group==1.5] = 2
AMBI$group[AMBI$group==4.5] = 4
AMBI$group[AMBI$group==6] = 5

## Read output of TITAN, mark pure and reliable taxa with p<0.05

# IndVal scores - represented also in heatmap
morphoInd = matrix(ncol = length(trainOn), nrow = dim(morpho)[2], 
                dimnames=list(names(morpho),trainOn))

# Groups for TITAN2 identified indicators, 
# reclassified as I -- V using quantile splines
morphoGrp = matrix(ncol = length(trainOn)+1, nrow = dim(morpho)[2], 
                dimnames=list(names(morpho),c(trainOn,"AMBIGroup")))

# Change points as determined by TITAN
morphoCP = matrix(ncol = length(trainOn), nrow = dim(morpho)[2], 
       dimnames=list(names(morpho),trainOn))

for (t in trainOn) {
  
  print(paste("---",t,"---"))
  
  included_samples = (!is.na(mds.nr.dna[,t]) & mds.nr.dna$Environment=="Estuarine")
  md.train = mds.nr.dna[included_samples,]
  row.names(md.train) = gsub("2017A","2017",row.names(md.train))
  row.names(md.train) = gsub("_PS","",row.names(md.train))
  table(row.names(md.train) %in% row.names(morpho))
  md.train = md.train[(row.names(md.train) %in% row.names(morpho)),]
  morpho.train = morpho[row.names(md.train),]

  ## Read TITAN2 output and select pure and reliable indicators
  oind = read.csv(paste("TITAN2/TITAN2_results/morpho_TITAN_all_",t,".csv", sep=""),header=T,row.names=1)
  oindS = oind[oind$reliability>=.95 & oind$purity>=.95,] # & oind$obsiv.prob<.05,]
  
  ## Use spline prediction
  mEGs = splinePredict(oindS,t, morpho.train, md.train, 
                         imageOutDir="img/sed/Splines/morpho/modelled")
  
  
  morphoGrp[,t] = mEGs$value # Does not match values!
  for (ti in c(1:dim(oindS)[1])){
    morphoInd[row.names(oindS)[ti],t] = oindS$IndVal[ti]
    morphoCP[row.names(oindS)[ti],t] = oindS$ienv.cp[ti]
  }
  
  ## Calculate BI values based on TITAN2 picked spline models
  
  indicators = mEGs[mEGs$value>0,,drop=F]
  if (dim(indicators)[1]==0) {
    print(paste("Error - no indicators for",t))
  } else {
    morpho_w_Ind = morpho.train[,row.names(indicators)]
  
    morpho_BI_values = rep(0,dim(morpho.train)[1])
    for (i in c(1:dim(indicators)[1])){
      w = getWeight(indicators[i,],bi=t)
      morpho_BI_values = morpho_BI_values + w*morpho_w_Ind[,i]
    }
    morpho_BI_values = morpho_BI_values/rowSums(morpho_w_Ind)
    
    ## Test correlation to real values (inside training set for now)
    
    preds = morpho_BI_values[!is.na(morpho_BI_values)]
    md.test = md.train[!is.na(morpho_BI_values),]
    
    pdf(paste("img/sed/Splines/morpho/",t,"_pred_morpho_noXVal",extra,".pdf",sep=""),width=8,height=5)
    plot_ml(data = preds, metadata = md.test, 
            xIndex = t, yIndex=t,
            aggreg = c("StationSeasonYear", "Estuary"), 
            title = "TITAN2+spline predicted BI values from morphotaxonomy")
    dev.off()
  }
}

# ----- Write output and make heatmaps -------

morphoGrp = as.data.frame(morphoGrp)
summary(morphoGrp)

## Make heatmap for taxon importance = IndVal

morphoInd[is.na(morphoInd)] = 0
totalmorphoInd = rowSums(morphoInd)
morphoInd_chart = morphoInd[order(totalmorphoInd,decreasing=TRUE),][c(1:50),]

pal <- colorRampPalette(c("black", "blue", "purple","red"))(n = 299)

pdf(paste("img/sed/TITAN2/Morphotaxon_importance_heatmap",extra,".pdf",sep=""),
    width=5.5,height=9)
heatmap(as.matrix(morphoInd_chart),col=pal,Colv = NA,
        Rowv = NA,revC=TRUE,scale="col",cexRow=.8,cexCol=.8,
        margins = c(5,10))
dev.off()

## Write TITAN2 summarised results
write.csv(morphoInd, "TITAN2/morpho_IndVal.csv")
write.csv(morphoGrp, "TITAN2/morpho_IndGroup.csv")
write.csv(morphoCP, "TITAN2/morpho_IndCP.csv")


