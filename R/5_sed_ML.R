## Anders Lanzén, 2019-09-22
## Trains RF classifier using molecular OTU and then taxon data, first using cross-validation
## then with all data, to determine importance. Heatmaps include EG assignments made using
## quantile regression splines. Methods and utility funcitons (sml_compo.R and ec_and_plot.R)
## are developed from scripts made available by Tristan Cordier (2019) 

## ****************************************************************************
## Note: Several operations in 1_sed_basic_analysis.R and 3_sed_TITAN.R 
## must have been run previous to this script!
## ****************************************************************************

source("R/utils/sml_compo.R")
source("R/utils/ec_and_plot.R")
require(irr)
require(ranger)

# --- RF training and cross-val, for each training variable -----

for (t in trainOn){
  
  included_samples = (!is.na(mds[,t]) & mds$Environment=="Estuarine" &
      mds$Nacid=="DNA")
  extra = "_DNA"

  md.train = mds[included_samples,]
  #md.train[md.train$AMBI==7,]$AMBI <- 6
  otus.train = otus.s.ra.f[included_samples,]
  table(row.names(md.train) == row.names(otus.train))
  
  taxa.train = taxa.ass.ra.f[included_samples,]
  
  otus.train = otus.train[,colSums(otus.train)>0]
  
  taxa.train = taxa.train[,colSums(taxa.train)>0]
  
  # ---- Predict RF classifier on all datasets with leave-one-out testing (on estuary name) ----
  
  preds <- sml_compo(otu_table = otus.train, metadata = md.train, index=t, 
                     cross_val = "Estuary", 
                     algo = "RF", optim_overfit = T)
  
  pdf(paste("img/sed/ML/OTUs_DNA/",t,"_RF_pred",extra,".pdf",sep=""),width=8,height=5)
  status <- plot_ml(data = preds, metadata = md.train, xIndex = t, yIndex=t,
                    aggreg = c("StationSeasonYear", "Estuary"), 
                    title = "All abundant OTUs (estuarine)")
  
  dev.off()
  
  
  # ----- Prediction using taxa instead of OTUs --------
  
  
  preds.taxa <- sml_compo(otu_table = taxa.train, metadata = md.train, index=t, cross_val = "Estuary",
                          algo = "RF", optim_overfit = T)

  pdf(paste("img/sed/ML/Taxa_DNA/",t,"_RF_pred_taxonomic",extra,".pdf",sep=""),width=8,height=5)
  status <- plot_ml(data = preds.taxa, metadata = md.train, xIndex = t, yIndex=t,
                    aggreg = c("StationSeasonYear", "Estuary"),
                    title = "All abundant taxa (estuarine)")
  dev.off()

  # Write predictions to disk
  predData = data.frame(row.names=row.names(md.train),
                        otu_preds=preds, taxa_preds=preds.taxa)
  write.table(predData,paste("MLPreds/",t,extra,".tsv",sep=""),
                sep="\t",quote=F, col.names=NA)
  
}


# ----- Most important OTUs for RF (uses all data to make new RF classifier)  --------
  # Selected according to Cordier et al (2018)

# Taxa importance across training sets matrices 

otuImp = matrix(ncol = length(trainOn), nrow = dim(otus.s.ra.f)[2], 
                dimnames=list(names(otus.s.ra.f),c(trainOn)))

taxaImp = matrix(ncol = length(trainOn), nrow = dim(taxa.ass.ra.f)[2],
                 dimnames=list(names(taxa.ass.ra.f),trainOn))


for (t in trainOn){
  
  # Select subset where data on training variable is available
  included_samples = (!is.na(mds[,t]) & mds$Environment=="Estuarine" & mds$Nacid=="DNA")
  extra = "_DNA" 
  
  md.train = mds[included_samples,]
  
  # uncomment to treat azoic samples as AMBI 6 effectively limiting to range 0--6 instead of 0--7
  #md.train[md.train$AMBI==7,]$AMBI <- 6
  
  otus.train = otus.s.ra.f[included_samples,]
  taxa.train = taxa.ass.ra.f[included_samples,]
  names(taxa.train) = make.names(names(taxa.train))
  
  tt = md.train[,t]
  
  ## For OTUs
  mod <- ranger(tt ~ ., data=otus.train, 
                mtry=floor(dim(otus.train)[2]/3), num.trees = 300, 
                importance= "impurity", write.forest = T)
  otuImp[,t] = mod$variable.importance
  
  # Most 30 important for graphs
  imp <- tail(sort(mod$variable.importance), 30)
  
  pdf(paste("img/sed/ML/OTUs_DNA/",t,"_OTU_importance",extra,".pdf",sep=""),width=5,height=6)
  p <- barplot(imp, horiz = T, xlab="Variable importance", 
               main=paste("OTU importance "), las=2, cex.names = 0.6)
  
  # get taxonomy
  tx <- otus.all[names(imp),"classification"]
  # get the last rank predicted
  bestTx = array(dim=length(tx))
  for (i in 1:length(tx)) bestTx[i] <- tail(unlist(strsplit(as.character(tx[i]), split=";", fixed=TRUE)), 1)
  # add to plot
  text(0, p, labels=bestTx, pos=4, cex=0.6)
  dev.off()

  ## For taxa
  modTax <- ranger(tt ~ ., data=taxa.train, 
                mtry=floor(dim(taxa.train)[2]/3), num.trees = 300, 
                importance= "impurity", write.forest = T)
  
  taxaImp[,t] = modTax$variable.importance
  
  imp <- tail(sort(modTax$variable.importance), 30)
  
  pdf(paste("img/sed/ML/Taxa_DNA/",t,"_Taxon_importance",extra,".pdf",sep=""),width=5,height=6)
  p <- barplot(imp, horiz = T, xlab="Variable importance", 
               main=paste("Taxon importance "), las=2, cex.names = 0.6)
  
  text(0, p, labels=names(imp), pos=4, cex=0.6)
  
  dev.off()

}

# --- Look at important OTUs and taxa across training sets ----

## Important OTUs

# Set NA values to zero
otuImp[is.na(otuImp)] = 0
otuImpNorm = decostand(t(otuImp),method="total")
totalOTUImp = colSums(otuImpNorm)
otuImp_chart = as.data.frame(otuImpNorm[,order(totalOTUImp,decreasing=TRUE)][,c(1:50)])
names(otuImp_chart)
rowSums(otuImp_chart)
#.16 (AMBI) -- .48 (PI HC) 
colSums(otuImp_chart)
# SWARM_1502  SWARM_2560   SWARM_286   SWARM_958   SWARM_989   SWARM_563   SWARM_125   SWARM_683 
# 0.41654595  0.21996020  0.13559174  0.11148085  0.11013606  0.09970018  0.09890934  0.09662508 
write.csv(otuImp_chart,"img/sed/ML/OTUs_DNA/otu_heatmap.csv")

# Add taxonomy information to most important OTUs
otuImp_tx = otus.all[names(otuImp_chart),"classification"]
bestImpTx=array(dim=50)
for (i in 1:50) bestImpTx[i] <- tail(unlist(strsplit(as.character(otuImp_tx[i]), split=";", fixed=TRUE)), 1)
bestImpTx[c(1:10)]
# [1] "Desulfobulbus"     "Sva1033"           "Woeseiaceae"       "Desulfobulbus"    
# [5] "Rhodospirillaceae" "Xanthomonadales"   "Woeseiaceae"       "Sva0725"          
# [9]  "Ca. Atribacteria" "Sva0071"  

## Important taxa
# Set NA values to zero
taxaImp[is.na(taxaImp)] = 0
taxaImpNorm  = decostand(as.data.frame(t(taxaImp)),method="total")
totalTaxaImp = colSums(taxaImpNorm)
taxaImp_chart = as.data.frame(taxaImpNorm[,order(totalTaxaImp,decreasing=TRUE)][,c(1:50)])
row.names(taxaImp_chart) = c("AMBI","PI","PI5y","PI Metals","PI HC","PI OM","Redox")
write.csv(taxaImp_chart,"img/sed/ML/Taxa_DNA/taxa_heatmap.csv")

rowSums(taxaImp_chart)
# AMBI        PI      PI5y PI Metals     PI HC     PI OM     Redox 
# 0.2538715 0.5509512 0.5812828 0.4250656 0.5397062 0.5218416 0.5357401 
colSums(taxaImp_chart[,c(1:10)])
# Mycobacterium       endosymbionts   Ca. Entotheonella               CA002 
# 0.43018994          0.21887113          0.19227449          0.16482760 
## ---- Make heatmaps -----

## Find ecogroups in spline predictions 

taxaImp_labels = taxaGrp[names(taxaImp_chart),]#c(8:10)]
for (i in c(1:50)){
  taxaImp_labels[i,] = as.character(as.roman(t(taxaImp_labels[i,])))
}

otuImp_labels = otuGrp[names(otuImp_chart),]#-8]
for (i in c(1:50)){
  otuImp_labels[i,] = as.character(as.roman(t(otuImp_labels[i,])))
  #otuImp_labels[i,is.na(otuImp_labels[i,])] = 0
  
}

MLOTUs = names(otuImp_chart)
names(otuImp_chart) = paste(gsub("SWARM_","",names(otuImp_chart)),bestImpTx)

#using a red and blue colour scheme without traces and scaling each row
library("RColorBrewer")
library("gplots")

pal <- colorRampPalette(c("black", "blue", "purple","red","yellow"))(n = 256)

pdf(paste("img/sed/ML/OTU_importance_heatmap",extra,".pdf",sep=""),
    width=6.2,height=10)
heatmap.2(t(sqrt(as.matrix(otuImp_chart))),col=pal,dendrogram="none",Colv = NA,
          Rowv = NA,revC=F,scale="none",cexRow=.7,cexCol=.8,
          margins = c(5,10),
          density.info="none",trace="none",
          cellnote = otuImp_labels,
          notecol="white",key=T,notecex=.6,
          keysize=2.2,
          key.par=list(mar=c(6,1,12,2)),
          
)
dev.off()


pdf(paste("img/sed/ML/Taxon_importance_heatmap",extra,".pdf",sep=""),
    width=6,height=10)
heatmap.2(t(sqrt(as.matrix(taxaImp_chart))),col=pal,dendrogram="none",Colv = NA,
          Rowv = NA,revC=F,scale="none",cexRow=.7,cexCol=.8,
          margins = c(5,10),
          density.info="none",trace="none",
          cellnote = taxaImp_labels,
          notecol="white",key=T,notecex=.6,
          keysize=2.2,
          key.par=list(mar=c(6,1,12,2)),
          
)
dev.off()


# # ----- Reference plots, lit. based microgAMBI etc. --------

for (t in trainOn){
  
  included_samples = (!is.na(mds[,t]) & mds$Environment=="Estuarine" & 
                        mds$Nacid=="DNA")
  extra = "_DNA"
  
  md.train = mds[included_samples,]
  
  
  if (t!="AMBI"){
    # v AMBI ref plot
    md.train.wAMBI = md.train[!is.na(md.train$AMBI),]
    pdf(paste("img/sed/BI_corrs/",t,"_v_AMBI",extra,".pdf",sep=""),height=5,width=8)
    plot_ml(data = md.train.wAMBI$AMBI, metadata = md.train.wAMBI, 
            xIndex = t, yIndex="AMBI",
            aggreg = c("StationSeasonYear", "Estuary"), 
            title = paste("AMBI v",t))
    dev.off()
  }
  
  if (t!="microgAMBI"){
    # v ugAMBI ref plot
    pdf(paste("img/sed/BI_corrs/",t,"_v_ugAMBI",extra,".pdf",sep=""),height=5,width=8)
    plot_ml(data = md.train$microgAMBI, metadata = md.train, 
            xIndex = t, yIndex="microgAMBI",
            aggreg = c("StationSeasonYear", "Estuary"), 
            title = paste("microgAMBI v",t))
    
    dev.off()
  }
}



