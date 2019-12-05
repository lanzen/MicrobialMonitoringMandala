setwd("~/projects/uGAMBI_RC/MBC_2019")

#library(labdsv)

# Tristan's stuffs 
#source("../R/computeIndexes.R")
source("R/sml_compo.R")
#source("../R/sml_compo_hyperparameter_fitting.R")
source("R/plot_ml.R")

require(irr)


## Data frames otus.s and mds must exist!

mds$StationSeasonYear = as.factor(paste(mds$Station,mds$Season, mds$Year,sep=""))
mds$col_plot = as.numeric(mds$Estuary)

# --------- Try different combinations of modelling ------

# trainOn = c("piMetals", "piHC")

for (t in trainOn){
  
  # included_samples = (!is.na(mds[,t]) & !is.na(mds$SalClass) & mds$Environment=="Estuarine")
  # extra=""
  # 
  included_samples = (!is.na(mds[,t]) & mds$Environment=="Estuarine" & mds$Nacid=="DNA")
  extra = "_DNA_only"
  # 
  # included_samples = (!is.na(mds[,t]) & mds$Environment!="Offshore")
  # extra = "_w_Coastal"
  # included_samples = (!is.na(mds[,t]))
  # extra = "_all_samples"
  
  md.train = mds[included_samples,]
  otus.train = otus.s.ra.f[included_samples,]
  
  taxa.train = taxa.ass.ra.f[included_samples,]
  
  otus.train = otus.train[,colSums(otus.train)>0]
  dim(otus.train) 
  
  # AMBI (est only): 149 samples, 4660 OTUs
  # AMBI (est, coast and offshore): 163 samples, 4685 OTUs
  # AMBI (est, coast and offshore): 177 samples, 4694 OTUs (performs v poorly)
  # AMBI (est): DNA only: 124 samples
  # PI (est): 188 samples, 4707 OTUs
  # PI (est) DNA only: 160 samples
  # PI 5/10y avg (est): 215 samples, 4707 (4) OTUs
  # PI avg (est), DNA only: 187 samples, 4702 OTUs
  # PI Metals DNA Only: 157 samples, 4701 OTUs
  # Redox (est): 184 sampels, 4706 OTUs
  # Redox (est), DNA only: 184 samples
  # Redox (est),  5 or 10y avg: 215 samples
  # OM (est): 188 samples
  # (OM (est) 5 or 10y avg: 217 samples)
  
  taxa.train = taxa.train[,colSums(taxa.train)>0]
  dim(taxa.train)
  # est only (regardless of index): 925 taxa
  
  # ------ Add more parameters as if they were OTUs -----
 
  # extra="_w_Shannon"
  # otus.train$shannon = div.s[included_samples,]$H #(almost same as richness)
  
  #otus.train = decostand(otus.train, method="hell") #(no effect whatsoever)
  # 
  #extra=extra+"_w_SalinityClass"
  # otus.train$Salinity = md.train$SalClass
  # 
  # extra="_w_depth"
  # otus.train$Depth = md.train$PointDepthBM
  
  # extra="_w_Mud"
  # otus.train$Mud = md.train$Mud

  # ---- Predict RF classifier on all datasets with leave-one-out testing (on estuary name) ----
  
  preds <- sml_compo(otu_table = otus.train, metadata = md.train, index=t, 
                     cross_val = "Estuary", 
                     algo = "RF", optim_overfit = T)
  
  pdf(paste("img/sed/ML/OTUs_DNA/",t,"_RF_pred",extra,".pdf",sep=""),width=8,height=5)
  status <- plot_ml(data = preds, metadata = md.train, index = t,
                    aggreg = c("StationSeasonYear", "Estuary"), 
                    title = "All abundant OTUs (estuarine)")
  
  dev.off()
  
  # ----- Prediction using taxa instead of OTUs --------
  
  
  preds.taxa <- sml_compo(otu_table = taxa.train, metadata = md.train, index=t, cross_val = "Estuary", 
                          algo = "RF", optim_overfit = T)
  
  pdf(paste("img/sed/ML/",t,"_RF_pred_taxonomic_rerun",extra,".pdf",sep=""),width=8,height=5)
  status <- plot_ml(data = preds.taxa, metadata = md.train, index = t,
                    aggreg = c("StationSeasonYear", "Estuary"), 
                    title = "All abundant taxa (estuarine)")
  
  dev.off()
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
  extra = "_DNA_only" 
  # included_samples = (!is.na(mds[,t]) & mds$Environment=="Estuarine")
  # extra=""
  # 
  md.train = mds[included_samples,]
  
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
  
  pdf(paste("img/sed/ML/",t,"_OTU_importance",extra,".pdf",sep=""),width=5,height=6)
  p <- barplot(imp, horiz = T, xlab="Variable importance", 
               main=paste("OTU importance "), las=2, cex.names = 0.6)
  
  # get taxonomy
  tx <- otus.t[names(imp),"classification"]
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
  
  pdf(paste("img/sed/ML/",t,"_Taxon_importance",extra,".pdf",sep=""),width=5,height=6)
  p <- barplot(imp, horiz = T, xlab="Variable importance", 
               main=paste("Taxon importance "), las=2, cex.names = 0.6)
  
  text(0, p, labels=names(imp), pos=4, cex=0.6)
  
  dev.off()

}

# --- Look at important OTUs and taxa across training sets ----

## Important OTUs

# Set NA values to zero
otuImp[is.na(otuImp)] = 0
otuImpNorm = decostand(as.data.frame(t(otuImp)),method="total")
totalOTUImp = colSums(otuImpNorm)
otuImp_chart = otuImpNorm[,order(totalOTUImp,decreasing=TRUE)][,c(1:50)]
rowSums(otuImp_chart)
# AMBI          PI     PI5yAvg       Redox  Redox5yAvg Redox10yAvg          OM     OM5yAvg    OM10yAvg 
# 0.3034950   0.4442061   0.5073967   0.4528319   0.6566396   0.6371903   0.4756302   0.6260022   0.6399512 
colSums(otuImp_chart[,c(1:10)])
# SWARM_238 SWARM_1225 SWARM_2094 SWARM_2412   SWARM_85 SWARM_1060  SWARM_374  SWARM_931 SWARM_1888 SWARM_1234 
# 0.39225923 0.30179769 0.23973078 0.17991454 0.16055340 0.13691709 0.13471769 0.12511186 0.12161234 0.09940919
otuImp_tx = otus.t[names(otuImp_chart),"classification"]
bestImpTx=array(dim=50)
for (i in 1:50) bestImpTx[i] <- tail(unlist(strsplit(as.character(otuImp_tx[i]), split=";", fixed=TRUE)), 1)
bestImpTx[c(1:10)]
# [1] "Woeseiaceae"       "Desulfobulbus"     "Rhodospirillaceae" "Draconibacterium"  "Woeseiaceae"      
# [6] "OM1 clade"         "Rhodospirillaceae" "Rhodospirillaceae" "Sva0725"           "Rhodospirillaceae"
# 4xRhodospirillaceae!

## Important taxa
# Set NA values to zero
taxaImp[is.na(taxaImp)] = 0
taxaImpNorm  = decostand(as.data.frame(t(taxaImp)),method="total")
totalTaxaImp = colSums(taxaImpNorm)
taxaImp_chart = taxaImpNorm[,order(totalTaxaImp,decreasing=TRUE)][,c(1:50)]
rowSums(taxaImp_chart)
# AMBI          PI     PI5yAvg       Redox  Redox5yAvg Redox10yAvg 
# 0.2159029   0.4962753   0.5678478   0.4761626   0.6239257   0.6394345 
# OM     OM5yAvg    OM10yAvg 
# 0.4400623   0.5664675   0.5899634 

colSums(taxaImp_chart[,c(1:10)])
# Sva0996 marine group                CA002     Acidimicrobiales    Ca. Entotheonella                 TK85 
# 0.4826956            0.3621616            0.2622194            0.2380116            0.1956823 
# Lweinellacea    Rhodospirillaceae        Mycobacterium         Trueperaceae      Marinilabiaceae 
# 0.1888974            0.1842415            0.1794176            0.1472758            0.1310021 

## Make heatmaps
names(otuImp_chart) = paste(gsub("SWARM_","",names(otuImp_chart)),bestImpTx)

#using a red and blue colour scheme without traces and scaling each row
library("RColorBrewer")
library("gplots")
pal <- colorRampPalette(c("black", "blue", "purple","red"))(n = 299)

pdf(paste("img/sed/ML/OTU_importance_heatmap",extra,".pdf",sep=""),
    width=5.5,height=9)
heatmap((t(as.matrix(otuImp_chart))),col=pal,Colv = NA,
            Rowv = NA,revC=TRUE,scale="col",cexRow=.8,cexCol=.8,
        margins = c(0,10))
#breaks=seq(0,.15,length.out = 300),trace="none")
dev.off()

pdf(paste("img/sed/ML/Taxon_importance_heatmap",extra,".pdf",sep=""),
    width=5.5,height=9)
heatmap((t(as.matrix(taxaImp_chart))),col=pal,Colv = NA,
        Rowv = NA,revC=TRUE,scale="col",cexRow=.8,cexCol=.8,
        margins = c(5,10))
dev.off()


# ----- Reference plots, lit. based microgAMBI etc. --------

pdf("img/sed/ML/microgAMBI_v_AMBI.pdf",height=5,width=8)
plot_ml(data = md.train$microgAMBI, metadata = md.train, index = t,
        aggreg = c("StationSeasonYear", "Estuary"), 
        title = "microgAMBI")

dev.off()


# ------------- SOMs ------------

# ---- SOM: Grid generation :----

library(mxnet)

grid_ra <- data.frame(array(NA, c(10,6)))
dimnames(grid_ra)[[2]] <- c("xdim", "ydim", "topo", "alpha", "radius", "xweight")

## generating random hyperparameters within intervals for SOM
for (i in 1:100)
{
  grid_ra[i,"xdim"] <- sample(1:10, 1)
  grid_ra[i,"ydim"] <- sample(1:10, 1)
  grid_ra[i,"topo"] <- c("rectangular", "hexagonal")[sample(1:2, 1)]
  grid_ra[i,"alpha"] <- paste("c(0.05,0.01)") 
  grid_ra[i,"radius"] <- runif(1, 0.2,2)
  grid_ra[i,"xweight"] <- runif(1, 0.2,1)
}
