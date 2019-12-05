## Anders Lanzén, 2019-09-22
## Compares RNA and DNA data from same samples. Generates several negative results not
## mentioned in the companion manuscript (Lanzén et al, 2020)

## ****************************************************************************************
## Several operations in 1_sed_basic_analysis.R must have been run previous to this script!
## ****************************************************************************************

require(vegan)
require(vioplot)
source('R/utils/diversity.r')
source('R/utils/correlationTests.r')
source('R/utils/taxaplot.R')

# -------- MicrogAMBI v AMBI, Variation between replicates, DNA v RNA -----

repCol = colour[mds$Rep=="A"]
lmMgA = lm(microgAMBI~AMBI,data=mds)
summary(lmMgA) # p<2E-16, adj R2=.35, 

pdf("img/sed/microgAMBI_vs_AMBI_simple.pdf",width=6,height = 6)
plot(mds$microgAMBI~mds$AMBI, col=colour, pch=as.numeric(mds$Environment),cex=.6,
     xlab="AMBI",ylab="microgAMBI")
for (i in c(1:dim(mds.repsA)[1])){
  lines(x=c(mds.repsA$AMBI[i],mds.repsB$AMBI[i]),
        y=c(mds.repsA$microgAMBI[i],mds.repsB$microgAMBI[i]),
        col=repCol[i])
}
abline(lmMgA,col="grey",lty=2)
legend("topleft",col=c(colU,"white",rep("black",3)), 
       legend=(c(as.vector(sort(unique(mdi$Estuary))),"","Estuary","Coast","Offshore")),
       pch=c(rep(1,15),2,1,3),cex=.6,ncol=1)

dev.off()

techRepDiff= abs(mds.repsA$microgAMBI - mds.repsB$microgAMBI)
summary(techRepDiff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.002645 0.009908 0.026800 0.025960 0.037060 0.053400 

lmMgD = lm(microgAMBI~AMBI,data=mds.nr.dna)
summary(lmMgD) #Adj R2=.41, mg =.84 + AMBI*.47 (almost identical to including all)
lmMgR = lm(microgAMBI~AMBI,data=mds.rnaComp)
summary(lmMgR) #Adj R2=-.002, p=.3, mg=2.6+AMBI*.11 (very poor)
lmMgDComp = lm(microgAMBI~AMBI,data=mds.dnaComp)
summary(lmMgDComp) #Adj R2=0.05, p=.06, mg=2.2+0.16*AMBI (almost but not sign.)

pdf("RNA_v_DNA/microgAMBI_vs_AMBI_RNA_vs_DNA_all.pdf",width=6,height = 6)
plot(mds$microgAMBI~mds$AMBI, pch=(15*(as.numeric(mds$Nacid)-1))+1,
     col=as.numeric(mds$Nacid)*2,cex=.8,
     xlab="AMBI",ylab="microgAMBI")
for (i in c(1:dim(mds.dnaComp)[1])){
  if (!is.na(mds.dnaComp$AMBI[i])){
    lines(x=c(mds.dnaComp$AMBI[i],mds.rnaComp$AMBI[i]),
          y=c(mds.dnaComp$microgAMBI[i],mds.rnaComp$microgAMBI[i]),col="purple")
  }
}
abline(lmMgA,col="grey",lty=2)
legend("topleft",pch=c(1,16),col=c(2,4),legend=c("DNA","RNA"))
dev.off()

# Is there a difference between how mgAMBI is predicted btw DNA v RNA?
RNAmgDiff= mds.rnaComp$microgAMBI - mds.dnaComp$microgAMBI
summary(RNAmgDiff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.66900 -0.03741  0.30530  0.19940  0.48570  1.92700 
wilcox.test(RNAmgDiff) #p=8e-4
# --> Yes, in average RNA predicts higher mgAMBI than DNA!

# Does RNA or DNA predicted mgAMBI better resemble AMBI?
RNA_DNA_mgAMBI_v_AMBI = data.frame(RNADiff = mds.rnaComp$microgAMBI - mds.rnaComp$AMBI,
                                   DNADiff = mds.dnaComp$microgAMBI - mds.dnaComp$AMBI,
                                   RNADiffAbs = abs(mds.rnaComp$microgAMBI - mds.rnaComp$AMBI),
                                   DNADiffAbs = abs(mds.dnaComp$microgAMBI - mds.dnaComp$AMBI))

summary(RNA_DNA_mgAMBI_v_AMBI)
# DNADiff lower than RNADiff but absolute difference about the same

wilcox.test(RNA_DNA_mgAMBI_v_AMBI$RNADiff,
            RNA_DNA_mgAMBI_v_AMBI$DNADiff,alternative="less")
wilcox.test(RNA_DNA_mgAMBI_v_AMBI$RNADiffAbs,
            RNA_DNA_mgAMBI_v_AMBI$DNADiffAbs,alternative="less") 

# --> Both underpredict but RNA less, so closer to AMBI, but far from significant (p=.6)

# ------- Alphadiv correlation, DNA v RNA -----

# Are ther any differences in correlations w regard to DNA v RNA diversity and paramters

printANOVA(mds.dnaComp[,c("Environment","Estuary","Season")], div.s.dnaComp[,-1], .05)
printANOVA(mds.rnaComp[,c("Environment","Estuary","Season")], div.s.rnaComp[,-1], .05)
# Su significantly less div. v winter for RNA (Simpson and Chao1) but not DNA (p=0.003)

printVS(div.s.dnaComp[,-1], mds.dnaComp[,vars], a=.01)
printVS(div.s.rnaComp[,-1], mds.rnaComp[,vars], a=.01)

#              DNA p   Adj.R2 |  RNA p   Adj R2
# RS~mgAMBI                      6E-7    .33
# RS~AMBI5y                      8E-4    .16  <---
# RS~Ni                          2E-6    .32
# RS~Redox5y                     .003    .12
# H~mgAMBI     1E-4   .21        5E-7    .34
# H~AMBI5y                       .008    .09
# H~Ni         .009   .10        3E-9    .45 <---
# H~Cr                           3E-4    .19
# H~Redox5y    .008   .10        8E-4    .16
# J'~mgAMBI    5E-5   .23        3E-7    .35 <---
# J'~PI5y      .002   .14
# J'~PIHC      .008   .10
# J'~Redox5y   9E-4   .16        4E-4    .18 <---
# J'~Ni        .007   .10        3E-9    .45 <---
# J'~Cr                          3E-5    .25 <---
# D~mgAMBI     6E-5   .23       7E-4     .34
# D~Ni        .009    .10        1E-8    .43
# D~Cr                           5E-8    .39
# D~Redox5y    .006   .11


# Test if microgAMBI is better described by RNA J' vs. DNA J'
mR = lm(mds.dnaComp$microgAMBI ~div.s.rnaComp$J)
mD = lm(mds.dnaComp$microgAMBI ~div.s.dnaComp$J)
anova(mR,mD) #RNA J' describes microgAMBI slightly better regardless of based on DNA or RNA
19
div.s.RNAvDNA = (div.s.rnaComp/div.s.dnaComp)[,-1]
summary(div.s.RNAvDNA)
# RS ranging from 0.64 to 1.2 and others similar

lmDRRS=lm(div.s.dnaComp$Rarefied.richness~div.s.rnaComp$Rarefied.richness)
summary(lmDRRS) #p=1E-9, Adj R2 .49
plot(div.s.dnaComp$Rarefied.richness~div.s.rnaComp$Rarefied.richness,
     xlim=c(2000,4000),ylim=c(2000,4000),xlab="RS (RNA)",ylab="RS (DNA)")
abline(1,1)

pdf("img/sed/Alphadiv/microgAMBI_vs_J_RNA_vs_DNA_all.pdf",width=6,height = 6)
plot(div.s$J ~ mds$microgAMBI, 
     xlab="microgAMBI",ylab="Pielou evenness",
     col=as.numeric(mds$Nacid)*2,cex=.8)
for (i in c(1:dim(mds.dnaComp)[1])){
  if (!is.na(mds.dnaComp$AMBI[i])){
    lines(x=c(mds.dnaComp$microgAMBI[i],mds.rnaComp$microgAMBI[i]),
          y=c(div.s.dnaComp$J[i],div.s.rnaComp$J[i]),col="purple")
  }
}
abline(lmMgA,col="grey",lty=2)

legend("bottomleft",ncol=2,pch=c(1,16),col=c(2,4),legend=c("DNA","RNA"))
dev.off()


pdf("img/sed/Alphadiv/PI_vs_J_RNA_vs_DNA_all.pdf",width=6,height = 6)
plot(div.s$J ~ mds$PI, 
     xlab="PI",ylab="Pielou evenness",
     col=as.numeric(mds$Nacid)*2,cex=.8)
for (i in c(1:dim(mds.dnaComp)[1])){
  if (!is.na(mds.dnaComp$AMBI[i])){
    lines(x=c(mds.dnaComp$PI[i],mds.rnaComp$PI[i]),
          y=c(div.s.dnaComp$J[i],div.s.rnaComp$J[i]),col="purple")
  }
}
abline(lmMgA,col="grey",lty=2)

legend("bottomleft",ncol=1,pch=c(1,16),col=c(2,4),legend=c("DNA","RNA"),cex=.8)
dev.off()

# ------- ANOSIM RNA v DNA -------------------
# Also see Giner et al (2019) where this comparison was significant in Malaspina data

mds.RNAvDNA = mds[mds$Comp=="D" | mds$Comp=="R",]
div.s.dnaComp = div.s[mds$Comp=="D" | mds$Comp=="R",]
anosim(otus.s.ra.f.RNAvDNA,mds.RNAvDNA$Comp)
# R=.003, p=0.3

# ------- adonis and Mantel RNA v DNA ---------

# Is there a difference in seasonal difference DNA v RNA?
anosim(otus.s.ra.f.rnaComp,mds.rnaComp$Season)
# R=.36, p=.003
anosim(otus.s.ra.f.dnaComp,mds.dnaComp$Season)
# R=.36, p=.002
# Not really

# Is there a difference between how DNA v RNA is structured with AMBI?
adonis(otus.s.ra.f.rnaComp[!is.na(mds.rnaComp$AMBI),]~AMBI,data=mds.rnaComp[!is.na(mds.rnaComp$AMBI),])
# NS     
adonis(otus.s.ra.f.dnaComp[!is.na(mds.dnaComp$AMBI),]~AMBI,data=mds.dnaComp[!is.na(mds.dnaComp$AMBI),])
# p=0.04m R2=.04 -> DNA OTU composition slightly better correlated with AMBI but not well

# Is there a difference between how DNA v RNA is structured with PI?
adonis(otus.s.ra.f.rnaComp[!is.na(mds.rnaComp$PI),]~PI,data=mds.rnaComp[!is.na(mds.rnaComp$PI),])
#R2=0.07***
adonis(otus.s.ra.f.dnaComp[!is.na(mds.dnaComp$PI),]~PI,data=mds.dnaComp[!is.na(mds.dnaComp$PI),])
#R2=.07***

adonis(otus.s.ra.f.rnaComp[!is.na(mds.rnaComp$Redox5yAvg),]~Redox5yAvg,data=mds.rnaComp[!is.na(mds.rnaComp$Redox5yAvg),])
#R2=.06***

adonis(otus.s.ra.f.dnaComp[!is.na(mds.dnaComp$Redox5yAvg),]~Redox5yAvg,data=mds.dnaComp[!is.na(mds.dnaComp$Redox5yAvg),])
#R2=.07***

mantel(vegdist(otus.s.ra.f.rnaComp), vegdist(otus.s.ra.f.dnaComp))
# Mantel statistic r: 0.745
# Significance: 0.001 



# ------- Taxa plots ------

taxa.dnaComp = taxa.s[,c(TRUE, TRUE, mds$Comp=="D")]
taxa.rnaComp = taxa.s[,c(TRUE, TRUE, mds$Comp=="R")]
table(row.names(mds.dnaComp) == names(taxa.dnaComp)[-c(1:2)])
table(row.names(mds.rnaComp) == names(taxa.rnaComp)[-c(1:2)])
table(mds[names(taxa.rnaComp),]$StationSeasonYear == mds[names(taxa.dnaComp),]$StationSeasonYear)

nonZero = (rowSums(taxa.rnaComp[,-c(1:2)])>0 | rowSums(taxa.dnaComp[,-c(1:2)])>0)
taxa.dnaComp = taxa.dnaComp[nonZero,]
taxa.rnaComp = taxa.rnaComp[nonZero,]

mds.comp = droplevels(rbind(mds.rnaComp,mds.dnaComp))
taxa.comp = cbind(taxa.rnaComp,taxa.dnaComp[,-c(1:2)])

grouping_info<-data.frame(row.names=row.names(mds.comp), mds.comp$Nacid)
ranks = data.frame(rank=c("domain","superkingdom","kingdom","phylum","class","order","family","genus","species"),
                   levels=c(3,11,26,26,26,26,26,26,26))

for (i in c(1:dim(ranks)[1])){
  level=toString(ranks$rank[i])
  levelTaxa = as.data.frame(t(taxa.comp[taxa.comp$Rank==level,-c(1:2)]))
  pdf(paste("RNA_v_DNA/img/",level,"_RNA_v_DNA.pdf",sep=""),height=7,width=)
  taxaplot(ranks$levels[i],grouping_info,levelTaxa)
  dev.off()
}


# ------- "Activity ratio" (rRNA/rDNA) ------

rRNARatio.taxa.log2 = log2((taxa.rnaComp[,-c(1:2)]) / (taxa.dnaComp[,-c(1:2)]))
plot(t(rRNARatio.taxa.log2["Mycobacterium",])~mds.dnaComp$PI)
plot(t(rRNARatio.taxa.log2["Ca. Nitrosopumilales",])~mds.dnaComp$PI)
plot(t(rRNARatio.taxa.log2["Syntrophaceae",])~mds.dnaComp$PI)

rRNARatio.otu.log2 = log2((otus.s.ra.f.rnaComp) / (otus.s.ra.f.dnaComp))
# rRNARatio.otu.log2[otus.s.ra.f.rnaComp == 0 & otus.s.ra.f.dnaComp > 0] = -10
# rRNARatio.otu.log2[otus.s.ra.f.dnaComp == 0 & otus.s.ra.f.rnaComp > 0] = 10
plot(rRNARatio.otu.log2[,"SWARM_1502"]~mds.dnaComp$piHC)
plot(rRNARatio.otu.log2[,"SWARM_125"]~mds.dnaComp$piOM)
plot(rRNARatio.otu.log2[,"SWARM_307"]~mds.dnaComp$AMBI)

table(row.names(taxa.ass.ra.f)==row.names(mds))

taxa.rnaComp.ass= taxa.ass.ra.f[mds$Comp=="R",]
taxa.dnaComp.ass = taxa.ass.ra.f[mds$Comp=="D",]

rRNARatio.ass.log2 = log2( taxa.rnaComp.ass / taxa.dnaComp.ass)

# Make heatmap with ML OTUs

otuActChart.ml = rRNARatio.otu.log2[,MLOTUs]
row.names(otuActChart.ml) = paste(mds.dnaComp$Station,mds.dnaComp$Season,mds.dnaComp$Year)
otuActChart.ml = otuActChart.ml[order(mds.dnaComp$PI),]

library("RColorBrewer")
library("gplots")
pal <- colorRampPalette(c("red", "yellow", "white","green", "blue"))(n = 299)

pdf("RNA_v_DNA/img/Ratio_heatmap_ML_OTUs.pdf", width=8,height=8)
heatmap(as.matrix(otuActChart.ml),Colv = NA,col=pal,
        Rowv = NA,revC=TRUE,scale="none",
        margins = c(8,8))
dev.off()

# Make heatmap with TITAN OTUs 

otuActChart.t = rRNARatio.otu.log2[,titanOTUs]
row.names(otuActChart.t) = paste(mds.dnaComp$Station,mds.dnaComp$Season,mds.dnaComp$Year)
otuActChart.t = otuActChart.t[order(mds.dnaComp$PI),]

pdf("RNA_v_DNA/img/Ratio_heatmap_TITAN_OTUs.pdf", width=8,height=8)
heatmap(as.matrix(otuActChart.t),Colv = NA,col=pal,
        Rowv = NA,revC=TRUE,scale="none",
        margins = c(8,8))
#breaks=seq(0,.15,length.out = 300),trace="none")
dev.off()


# Make heatmap with ML taxa

taxActChart.ml = rRNARatio.ass.log2[,names(taxaImp_chart)]
row.names(taxActChart.ml) = paste(mds.dnaComp$Station,mds.dnaComp$Season,mds.dnaComp$Year)
taxActChart.ml = taxActChart.ml[order(mds.dnaComp$PI),]

pdf("RNA_v_DNA/img/Ratio_heatmap_ML_taxa.pdf", width=8,height=8)
heatmap(as.matrix(taxActChart.ml),Colv = NA,col=pal,
        Rowv = NA,revC=TRUE,scale="none",
        margins = c(10,7))
#breaks=seq(0,.15,length.out = 300),trace="none")
dev.off()


# Make heatmap with Spline taxa

taxActChart.t = rRNARatio.ass.log2[,titanTaxa]
row.names(taxActChart.t) = paste(mds.dnaComp$Station,mds.dnaComp$Season,mds.dnaComp$Year)
taxActChart.t = taxActChart.t[order(mds.dnaComp$PI),]

pdf("RNA_v_DNA/img/Ratio_heatmap_TITAN_taxa.pdf", width=8,height=8)
heatmap(as.matrix(taxActChart.t),Colv = NA,col=pal,
        Rowv = NA,revC=TRUE,scale="none",
        margins = c(10,7))
#breaks=seq(0,.15,length.out = 300),trace="none")
dev.off()

# Make heatmap with top 50 orders and OTUs

rank="order"
rankTaxa = taxa.comp[taxa.comp$Rank==rank,]
toptaxa = row.names(rankTaxa)[order(rowSums(rankTaxa[,-c(1:2)]),decreasing = T)]
taxActChart = rRNARatio.taxa.log2[toptaxa[c(1:50)],]
names(taxActChart) = paste(mds.dnaComp$Station,mds.dnaComp$Season,mds.dnaComp$Year)
taxActChart = taxActChart[,order(mds.dnaComp$PI)]

pdf("RNA_v_DNA/img/Ratio_heatmap_top50Orders.pdf",width=8,height=8)
heatmap((t(as.matrix(taxActChart))),Colv = NA,col=pal,
        Rowv = NA,revC=TRUE,scale="none",
        margins = c(11,8))
dev.off()


# -------- RNA v DNA plot -----------

RNAvDNA = data.frame(rna = colMeans(taxa.rnaComp.ass)+1E-6, 
                     dna = colMeans(taxa.dnaComp.ass)+1E-6,
                     ratio = colMeans(taxa.rnaComp.ass) / colMeans(taxa.dnaComp.ass))

asMeanML = RNAvDNA[names(taxaImp_chart),]
asMeanSpline = RNAvDNA[names(taxaInd_chart),]

plot(rna~dna, data = RNAvDNA, log="xy", cex=.8,col="grey", xlab="DNA abundance",
     ylab="RNA abundance",xlim=c(1e-6,.2),ylim=c(1e-6,.2))
abline(0,1,lty=2)
legend("bottomright",box.lwd=0,legend=c("Top 50 ML","Top 50 IndVal"),
       pch=c(17,18),cex=.8)
points(rna~dna, asMeanSpline, pch=18,cex=1.2)
points(rna~dna, asMeanML, pch=17)

topAb = row.names(RNAvDNA)[order(RNAvDNA$dna,decreasing=T)]
topAb = topAb[topAb %in% row.names(asMeanML) | topAb %in% row.names(asMeanSpline)]
lowRatio = row.names(RNAvDNA)[order(RNAvDNA$ratio)]
lowRatio = lowRatio[lowRatio %in% row.names(asMeanML) | lowRatio %in% row.names(asMeanSpline)]


text(rna~dna, RNAvDNA[topAb[c(1:5)],], labels= topAb[c(1:5)],cex=.5, pos=4)
text(rna~dna, RNAvDNA[lowRatio[c(1:5)],], labels= lowRatio[c(1:5)],cex=.5, pos=4,col="red")

RNAvDNA[topAB[c(1:5)],]

RNAvDNA[lowRatio[c(1:5)],]

# rna          dna      ratio
# Christensenellaceae  4.480297e-05 0.0007203303 0.06089409
# NB1-n                5.024207e-04 0.0037559232 0.13353688
# Parasphingopyxis     2.478680e-05 0.0001621625 0.14759512
# FD035                9.216290e-05 0.0006155492 0.14834109
# Ca. Nitrosopumilales 6.765039e-04 0.0037296438 0.18116610

## Plot for OTUs

RNAvDNA.otus  = data.frame(rna = colMeans(otus.s.ra.f.rnaComp)+1E-6, 
                     dna = colMeans(otus.s.ra.f.dnaComp)+1E-6,
                     ratio = colMeans(otus.s.ra.f.rnaComp) / colMeans(otus.s.ra.f.dnaComp))

otuRatio = RNAvDNA.otus$ratio[is.finite(RNAvDNA.otus$ratio) ]
summary(otuRatio)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.0000   0.3766   0.8240   1.8440   1.7570 120.6000 
wilcox.test(otuRatio-1) 
#p=.013

asMeanML.otus = RNAvDNA.otus[MLOTUs,]
asMeanSpline.otus = RNAvDNA.otus[titanOTUs,]

indiRatio.otus = c(asMeanML.otus$ratio, asMeanSpline.otus$ratio)

summary(indiRatio.otus)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.3688  0.6822  0.9779  1.1430 13.2600 

wilcox.test(indiRatio.otus,otuRatio,mu=0, alternative="less") 
# p=.02


plot(rna~dna, data = RNAvDNA.otus, log="xy", cex=.6,col="grey", xlab="DNA abundance",
     ylab="RNA abundance",xlim=c(1e-6,.05),ylim=c(1e-6,.05))
abline(0,1,lty=2)
legend("bottomright",box.lwd=0,legend=c("Top 50 ML","Top 50 IndVal"),
       pch=c(17,18),cex=.8)
points(rna~dna, asMeanSpline.otus, pch=18,cex=1.2)
points(rna~dna, asMeanML.otus, pch=17)

topAb.otus = row.names(RNAvDNA.otus)[order(RNAvDNA.otus$dna,decreasing=T)]
topAb.otus = topAb.otus[topAb.otus %in% row.names(asMeanML.otus) | topAb.otus %in% row.names(asMeanSpline.otus)]
lowRatio.otus = row.names(RNAvDNA.otus)[order(RNAvDNA.otus$ratio)]
lowRatio.otus = lowRatio.otus[lowRatio.otus %in% row.names(asMeanML.otus) | lowRatio.otus %in% row.names(asMeanSpline.otus)]

tx <- otus.all[topAb.otus,"classification"]
topTx = array(dim=length(tx))
for (i in 1:length(tx)) topTx[i] <- tail(unlist(strsplit(as.character(tx[i]), split=";", fixed=TRUE)), 1)
topTx = paste(gsub("SWARM_","",topAb.otus) ,topTx)

tx <- otus.all[lowRatio.otus,"classification"]
lowTx = array(dim=length(tx))
for (i in 1:length(tx)) lowTx[i] <- tail(unlist(strsplit(as.character(tx[i]), split=";", fixed=TRUE)), 1)
lowTx = paste(gsub("SWARM_","",lowRatio.otus) ,lowTx)

text(rna~dna, RNAvDNA.otus[topAb.otus[c(1:5)],], labels= topTx[c(1:5)],cex=.6, pos=4)
text(rna~dna, RNAvDNA.otus[lowRatio.otus[c(1:5)],], labels= lowTx[c(1:5)],cex=.6, pos=4,col="red")

RNAvDNA.otus[lowRatio.otus[c(1:5)],]

# rna          dna      ratio
# SWARM_3149 1.000000e-06 2.646806e-05 0.00000000
# SWARM_5036 2.320101e-06 6.911979e-05 0.01937911
# SWARM_285  6.526108e-05 2.935499e-03 0.02189849
# SWARM_882  3.419454e-05 4.779885e-04 0.06959190
# SWARM_1410 7.925195e-06 7.032916e-05 0.09988864

rdo = RNAvDNA[order(RNAvDNA$ratio, decreasing = T),]
rdo[rdo$dna>1E-3,][c(1:10),]

# ------- NMDS RNA only -------

nmds.RNA = metaMDS(otus.s.ra.f.rnaComp)
ef.rna = envfit(nmds.RNA, mds.rnaComp[,c(16,21:30,33:34,39:44,50,54:55,73:75)],na.rm=T)
ef2.rna = envfit(nmds.RNA, mds.rnaComp[,c(17,19)],na.rm=T)

rnacol = rb[as.numeric(mds.rnaComp$Estuary)] 
pch.rna=as.numeric(mds$Environment)-1
pch.rna[mds$Season=="Su" & mds$Nacid=="DNA" & mds$Environment=="Estuarine"] = 5
pch.rna[mds$Season=="Su" & mds$Nacid=="RNA"] = 18
pch.rna[mds$Season=="Wi" & mds$Nacid=="RNA" & mds$Environment=="Estuarine"] = 16
pch.rna[mds$Season=="Wi" & mds$Nacid=="RNA" & mds$Environment=="Coastal"] = 15

ordiplot(nmds.s,type="none",ylim=c(2,-2))
points(nmds.RNA,display="sites",col=rnacol,cex=.9,lwd=1,pch=pch.rna)
#text(nmds.s,display="sites",col=nmdscol, cex=0.2, labels=row.names(mds), pos=1)
plot(ef.rna,p.max=0.001,col="purple",cex=.7)
legend("bottomright",box.lwd=0.5, pch=c(1,5,0,2,3),
       legend=c("Estuarine Wi","Estuarine Su","Coastal","Offshore","Port"),
       ncol=1,cex=.6)
legend("bottomleft",box.lwd=0.5, legend=sort(unique(mds$Estuary)),col=colU,
       cex=.6,ncol=3,pch=1)

