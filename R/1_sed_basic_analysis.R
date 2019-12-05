## Anders Lanzén, 2019-09-22
## Reads OTUs adn metadata and calculates basic statistics

setwd("~/projects/uGAMBI_RC/Mandala_Sediments_Article1") #<- change to correct local directory
otudir = "CREST_Final_Prok"

require(vioplot)
require(vegan)
source('R/utils/filtering.R')
source('R/utils/diversity.r')
source('R/utils/correlationTests.r')
source('R/utils/taxaplot.R')
source("R/utils/calculatePIs.R")

# ----- Read OTU table and metadata ---------

## Choose variables to train on for ML and TITAN/Splines:
trainOn = c("AMBI", "PI", "PI5yAvg", "piMetals", "piHC", "piOM", "Redox")

## Read metadata
mds = read.table(file="metadata/metadata_final_S.tsv",header=T,row.names=1,sep="\t")
mds = mds[sort(row.names(mds)),]

mds$intertidal = mds$PointDepthBM==0
mds$Date = as.Date(mds$Date)

## Read OTU-table
otu.file = paste(otudir, "SWARM_table_final.tsv.gz", sep="/")
otus.all = read.delim(otu.file,row.names=1,header=T,sep="\t")
otus.all = otus.all[,order(names(otus.all))]

taxonomy=otus.all$classification
otus.s = as.data.frame(t(otus.all[,-1]))
otus.s = otus.s[,colSums(otus.s)>0]
taxonomy.s = taxonomy[colSums(otus.s)>0]

## Check agreement to metadata
table(row.names(otus.s) == row.names(mds))

# Read diversity results and check number of samples, reads and OTUs
dim(otus.s) # 316 samples, 185340 OTUs
sum(otus.s) # 10,416,505 reads

divTemp = read.csv("Diversity_stats/Diversity_classified_S.csv", row.names=1, header=T)
div.s = divTemp[,c(1,5:7,9:10)]

# Check agreement to otus
table(row.names(div.s) == row.names(mds))

summary(div.s)
# 10-138t reads, Rarefied richness = 1139 -- 4800, H' = 4.8 -- 8.3, J' = 0.61 -- 0.93

## ---- Metadata operations ------

mds = calculatePIs(mds)

# Parameters of factor type
mdf = mds[,c("Nacid","Run","Primer","Environment","Estuary","Station","Season",
             "Year","SalClass", "intertidal")]

# All continous parameters
vars = -c(1:14,45:49,51:52,60:72)
mdc = mds[,vars]

## ----- Groups of continuous variables ----

mdse = mds[mds$Environment=="Estuarine",]

# Continous parameters available for *all* samples (for envit / NMDS and similar)
md1 = mds[,c("Year","microgAMBI")]
md1e = mdse[,c("Year","microgAMBI", "PI","TotalPCB","TotalPAH","TotalDDT","TotalHCH","piHC",
               "OM","Gravel","Sand","Mud","piOM")]

# Parameters grp 2. (all but port and heavily contaminated, n=8 missing)
md2 = mds[,c("Salinity","SalClass","OM5yAvg","piOM5yAvg",
              "MaxDepthPM","PointDepthBM")]
md2e = mdse[,c("Salinity","SalClass","OM5yAvg","piOM5yAvg",
             "MaxDepthPM","PointDepthBM")] # 4 missing

# Parameters grp 3 (n=65 missing)
md3 = mds[,c("AMBI","AMBI")]
md3e = mdse[,c("AMBI","AMBI")] #61 missing

# Parameters grp 6 (n=25 missing)
md6= mds[,c("Redox","Redox")] 
md6e= mdse[,c("Redox","Redox")] # 9 missing

# Parameters grp 7 (n=20 missing)
md7=mds[,c("Cu","Pb","Ni","Cr","Zn","Hg","Cd","piMetals")] 
md7e=mdse[,c("Cu","Pb","Ni","Cr","Zn","Hg","Cd","piMetals")]  # 4 missing

# Parameters grp 8 (n=39 missing)
md8= mds[,c("Mn","Fe")]
md8e= mdse[,c("Mn","Fe")] # 22 missing

# Parameters grp 9 (n=4 missing)
md9= mds[,c("OM","Gravel","Sand","Mud","piOM")]


# Parameters grp 11 (n=16 missing)
md11= mds[,c("PI","PI5yAvg","TotalPCB","TotalPAH","TotalDDT","TotalHCH","piHC")]

# Parameters grp 12 (n=22 missing)
md12= mds[,c("WaterO2","WaterO2")]
md12e= mdse[,c("WaterO2","WaterO2")] # 6 missing

# Parameters grp 13 (n=42 missing)
md13= mds[,c("Q1SiteSal","Q2SiteSal","Q3SiteSal")]

# ------ Calculate relative abundance and filter rare OTUs ------
otus.s.ra.all = decostand(otus.s, method="total")
otus.s.ra.f1 = dropRareByMaxAbundance(otus.s.ra.all, 10/min(rowSums(otus.s)))
otus.s.ra.f = decostand(otus.s.ra.f1, method="total")
dim(otus.s.ra.f) # 5203 otus of >180k retained

taxonomy.s.ra.f = as.character(otus.all[names(otus.s.ra.f),"classification"])

otus.s.ra.bc = vegdist(otus.s.ra.f)

# Graph colours
rb = rainbow(17)
colour = rb[as.numeric(mds$Estuary)]
colU = rb[sort(as.numeric(unique(mds$Estuary)))]
colour[mds$Estuary == "None"] = "black"
colU[sort(unique(mds$Estuary)) == "None"] = "black"

# -------- Make subsets for removing replicates, DNA or RNA only ----------

# Subset without replicated samples for parameter correlation (removing B replicates, n=312)
mds.nr = mds[mds$Rep!="B",]
div.s.nr = div.s[mds$Rep!="B",]
otus.s.ra.f.nr = otus.s.ra.f[mds$Rep!="B",]

# Subset without replicated sample and with only DNA (for statsitics, removing RNA)

mds.nr.dna = mds.nr[mds.nr$Nacid=="DNA",]
mdf.nr.dna = mds.nr.dna[,c("Nacid","Run","Primer","Environment","Estuary","Station","Season",
             "Year","SalClass", "intertidal")]
div.s.nr.dna = div.s.nr[mds.nr$Nacid=="DNA",]
otus.s.ra.f.nr.dna = otus.s.ra.f.nr[mds.nr$Nacid=="DNA",]

# Subset of samples that have replicates (n=5)

mds.repsA = mds[mds$Rep=="A",]
mds.repsB = mds[mds$Rep=="B",]
div.s.repsA = div.s[mds$Rep=="A",]
div.s.repsB = div.s[mds$Rep=="B",]
otus.s.ra.f.repsA = otus.s.ra.f[mds$Rep=="A",]
otus.s.ra.f.repsB = otus.s.ra.f[mds$Rep=="B",]
otus.s.ra.f.repsBoth = rbind(otus.s.ra.f.repsA,otus.s.ra.f.repsB)

# Bray-Curtis dissimilarity between technical replicates from same sample
bcRepDiff = diag(as.matrix(vegdist(otus.s.ra.f.repsBoth))[c(1:5),c(6:10)])
summary(bcRepDiff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.111  0.156  0.157    0.155  0.172   0.177 

# DNA vs RNA subset for samples were both were sequenced, n=60
# All final RNA libs have DNA companion

mds.dnaComp = mds[mds$Comp=="D",]
mds.rnaComp = mds[mds$Comp=="R",]
div.s.dnaComp = div.s[mds$Comp=="D",]
div.s.rnaComp = div.s[mds$Comp=="R",]
otus.s.ra.f.dnaComp = otus.s.ra.f[mds$Comp=="D",]
otus.s.ra.f.rnaComp = otus.s.ra.f[mds$Comp=="R",]
otus.s.ra.f.RNAvDNA = rbind(otus.s.ra.f.dnaComp,otus.s.ra.f.rnaComp)

# Bray-Curtis for DNA v RNA from same sample
bcRNAvDNA_Diff = diag(as.matrix(vegdist(otus.s.ra.f.RNAvDNA))[c(1:60),c(61:120)])
write.table(as.matrix(vegdist(otus.s.ra.f.RNAvDNA))[c(1:60),c(61:120)],
            "RNA_v_DNA/BCDissimilarity.tsv", sep="\t",quote=F, col.names=NA)
summary(bcRNAvDNA_Diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.266  0.358  0.405   0.446  0.457    0.929 <-- EB05_S_Wi2013RNA
# This outlier is strange but likely results from RNA and DNA extracted and amplified different years


interSampleBC = as.numeric(unlist(vegdist(otus.s.ra.f.nr.dna)))
summary(interSampleBC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0118  0.681  0.802  0.777  0.892  0.996

wilcox.test(bcRepDiff,interSampleBC,alternative="less") # p=5E-5, replicates are more similar
wilcox.test(bcRNAvDNA_Diff,interSampleBC,alternative="less") # p<2E-16, DNA-RNA are more similar
wilcox.test(bcRNAvDNA_Diff,bcRepDiff) # p=2-4E-4 DNA-RNA are more different than tech reps 

grpinfo = data.frame(row.names=row.names(otus.s.ra.f.RNAvDNA), c(rep("DNA",60),rep("RNA",60)))
pdf("RNA_v_DNA/OTU_plot.pdf",width=20,height=12)
taxaplot(50,grpinfo,otus.s.ra.f.RNAvDNA)
dev.off()


# Coastal and Offshore only
oc = otus.s.ra.f.nr.dna[mds.nr.dna$Environment=="Coast" | mds.nr.dna$Environment=="Offshore",]
ocInter = as.numeric(unlist(vegdist(oc)))
summary(ocInter)


pdf("img/sed/BCDifs_viol_replication_Comp.pdf",width=5.7,height=5.7)

vioplot(bcRepDiff,bcRNAvDNA_Diff,interSampleBC,ocInter,
        names=c("Replicates","RNA v DNA","Intersample", "Littoral"))
dev.off()

# -------- Parameter correlations  ---------  

# Compare factors to continuous variables, in unreplicated DNA subset, Bonferroni correction
printANOVA(droplevels(mdf.nr.dna[,-c(1:3,5:6)]),
           mds.nr.dna[,vars],a=.05/(47*7))

# Compare continuous variables to each other, in unreplicated DNA subset, Bonferroni 
printAllvsAll(mds.nr.dna[,vars],a=.05/2850)
# mgAMBI ~ Redox: p=7E-16, R2=.25 (better than AMBI)
# AMBI ~ Redox: p=5E-9, R2=.17
# mgAMBI ~ Redox.5y: p<2E-16, R2=.47
# AMBI ~ Redox.5y: p=7E-14, R2=.26
# mgAMBI ~ OM: p=2E-9, R2=.13 (comparable to AMBI) (PI_OM: 4E-16, R2=.23)
# AMBI ~ OM: p=7E-7, R2=.12 (PI_OM: 2E-8, R2=.15)
# mgAMBI ~ OM.5y: p=9E-12, R2=.17 (comparable to AMBI, better than yearly)
# AMBI ~ OM.5y: p=6E-11, R2=.20 (PI_OM5y: 2E-8, R2=.15)
# mgAMBI ~ piOM: R2=.23 p=4E-16 (much better and more linear than OM)
# mgAMBI ~ piOM5y: R2=.23 p=8E-16 (better still)

# mgABMI ~ mud: p=1E-12, R2=0.18 
# AMBI ~ Mud: p=6E-9, R2=.16

# mgAMBI ~ PI: 1=7E-12, R2=0.19 (better than AMBI but 5y AMBI average best)
# AMBI ~ PI: p=2E-8, R2=.16
# mgAMBI ~ PI.5y: p=2E-16, R2=.27 (better than yearly / integrative, much better than AMBI)
# AMBI ~ PI.5y: p=1E-8, R2=.16

# mgAMBI ~ Ni: p=2E-8, R2=.12
# AMBI ~ Ni: p=2E-4, R2=.07

# mgAMBI ~ PI_Metals: p=3E-8, R2=.12
# AMBIy ~ PI_Metals: 7E-7, R2=.12

# mgAMBI ~ PI_HC: p=3E-8 R2=.12

# mgAMBI ~ salinity (neg): R2=.12, p=2E-8 (stronger for AMBI, R2=.14 same p)

# Very interesting that mgAMBI appears to correlate stronger with 5y averages for PI, OM 
# and redox. Though AMBI does also, except perhaps weaker and not(comparable to AMBI) for PI.

# Redox ~ OM (neg): p<2E-16, R2=.38 (less for redox 5 or 10y) (better fit with PI_OM:R2=.52)
# Redox ~ Mud (neg): p<2E-16, R2=.44 (5y Redox .40)
# Redox ~ PON (neg): p=7E-7, R2=.50
# Redox ~ TotalP (neg): p=4E-5, R2=.02, Bonferroni correction
# Redox ~ PI (neg):  p<2E-16, R2=.34 (5y 0.36) 
# Redox ~ piOM: p<2E-16 R2=.51 (very strong correlation, expected)
# Redox ~ piMetals: p=1E-11, R2=.18
# Redox ~ piHC: p<2E-16, R2=.37

# Mud ~ OM: p<2E-16, R2=.28 (much stronger, R2=.44 with PI_OM)
# PI ~ OM: p<2E-16, R2=.26 (5y: R2=0.23)
# OM ~ pIMetals: 1E-7, 0.12
# OM ~ pIHC: 9E-15, 0.22
# Mud ~ PI: p=1E-13, R2=.21 (5y similar R2=.21)
# Metals, PCBs etc to each other
# PI ~ Cd: p<2E-16, R2=0.62 (5y: R2=0.47) <--
# PI ~ Cu: p<2E-16, R2=0.53 (5y: R2=0.36) <--
# PI ~ Zn: p<2E-16, R2=0.55 (5y: R2=0.41) <--
# Other metals like Pb also well correlated (R2=0.2 -- 0.4)
# TotalPCB ~ PI: p<2E-16, R2=0.56 (PI5y 0.42) <-

# ------- Alphadiv correlation -----

# Compare alpha-diversity measures with continous physicochemical ones
printANOVA(droplevels(mdf.nr.dna), div.s.nr.dna, .05)

# Higher Chao1 with new primers but likely selection bias
# Coastal has lower RS v Estuarine (p=0.016)
# Much clearer for Chao1 (p=7E-7) and offshore too (3E-7)
# Su lower in RS than Wi (p=0.02), H' (p=0.04) and Chao1 (p=.004)
# RS and Chao1 follow laying S with salinity class gradient (1 highest, 2.5 local max, 4 lowest)
# Intertidal more diverse. Probably Nervion

de = div.s.nr.dna[mdf.nr.dna$Environment=="Estuarine",]
dm = div.s.nr.dna[mdf.nr.dna$Environment!="Estuarine",]

wilcox.test(de$Rarefied.richness,dm$Rarefied.richness)
# p= 0.02

me = mds.nr.dna[mds.nr.dna$Environment=="Estuarine",]
mm = mds.nr.dna[mds.nr.dna$Environment!="Estuarine",]
wilcox.test(me$AMBI,mm$AMBI) # p=1E-12
wilcox.test(me$microgAMBI,mm$microgAMBI) # p=8E-12
wilcox.test(me$Redox,mm$Redox) # p=2E-6


printVS(div.s.nr.dna, mds.nr.dna[,vars], a=.05/49)
#             p     Adj R2
# RS~Cu       6E-5    .06 (+other metals)
# RS~piMetal  1E-5    .08....
# RS~PI       2E-4    .06
# RS~PI5y     1E-4    .06
# H~PIMetals  7E-6    .08
# H'~PIHC     3E-4    .05
# H'~PI       9E-6    .08
# H'~PI5y     3E-6    .09
# H'~Cu       1E-4    .06
# J'~Cu       1E-4    .06
# J~PI        5E-6    .08
# J~PI5y      3E-7    .10
# J~PIMetals  2E-5    .07
# J~PIHC      8E-5    .06
# J~PI5y      2E-4    .06

summary(lm(div.s.nr.dna$H~mds.nr.dna$AMBI))
# R2 = .07, p=2E-4
summary(lm(div.s.nr.dna$J~mds.nr.dna$AMBI))
# R2 = .07, p=7E-5

###
repCol = colour[mds$Rep=="A"]


m1 = glm(div.s$J~mds$PI)

pdf("img/sed/Alphadiv/PI_vs_evenness.pdf",width=6,height = 6)
plot(div.s$J ~ mds$PI, col=colour, pch=as.numeric(mds$Environment),cex=.6,
     xlab="PI",ylab="Pielou evenness")
for (i in c(1:dim(mds.repsA)[1])){
  lines(x=c(mds.repsA$PI[i],mds.repsB$PI[i]),
        y=c(div.s.repsA$J[i],div.s.repsB$J[i]),
        col=repCol[i])
}
abline(m1,col="grey",lty=2)
dev.off()


# ------- ANOSIM and adonis  --------

anosim(otus.s.ra.f.nr.dna, grouping=mds.nr.dna$Station)
# R=.71 ***

anosim(otus.s.ra.f.nr.dna, grouping=mds.nr.dna$Estuary)
# R=.26***

anosim(otus.s.ra.f.nr.dna, grouping=mds.nr.dna$Environment)
# R=.39 ***

anosim(otus.s.ra.f.nr.dna, grouping=mds.nr.dna$Season)
# NS

anosim(vegdist(otus.s.ra.f.nr.dna[!is.na(mds.nr.dna$SalClass),]), 
       grouping=mds.nr.dna$SalClass[!is.na(mds.nr.dna$SalClass)])
#R=.21***

adonis(otus.s.ra.f.nr.dna~Environment+Season+Station+Year,data=mds.nr.dna)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Environment   3     8.519 2.83980 19.6095 0.10937  0.001 ***
# Season        1     2.125 2.12542 14.6766 0.02729  0.001 ***
# Station      64    40.472 0.63237  4.3667 0.51959  0.001 ***
# Year          1     0.563 0.56342  3.8905 0.00723  0.001 ***
# Residuals   181    26.212 0.14482         0.33652           
# Total       250    77.892                 1.00000     


# ----- NMDS -----

nmds.s = metaMDS(otus.s.ra.f)
#Run 20 stress 0.1377538 

# Correlations to sample positions in NMDS space modelled with envfit individually
# due to missing data differences
efs = list(envfit(nmds.s,md1),envfit(nmds.s,md2,na.rm=T),
        envfit(nmds.s,md3,na.rm=T),envfit(nmds.s,md5,na.rm=T),
        envfit(nmds.s,md6,na.rm=T),envfit(nmds.s,md7,na.rm=T),envfit(nmds.s,md8,na.rm=T),
        envfit(nmds.s,md9,na.rm=T),envfit(nmds.s,md11,na.rm=T),
        envfit(nmds.s,md12,na.rm=T),envfit(nmds.s,md13,na.rm=T))

svg("img/sed/NMDS_w_envfit.svg")
  
  ordiplot(nmds.s,type="none")
  for (i in c(1:3,5:8, 10:11)){
    plot(efs[[i]],col="purple",p.max=.001,cex=.5)
  }
  nmdscol=colour
  #nmdscol=as.numeric(mds$Nacid)*2
  pch=as.numeric(mds$Environment)-1
  pch[mds$Season=="Su" & mds$Nacid=="DNA" & mds$Environment=="Estuarine"] = 5
  pch[mds$Season=="Su" & mds$Nacid=="RNA"] = 18
  pch[mds$Season=="Wi" & mds$Nacid=="RNA" & mds$Environment=="Estuarine"] = 16
  pch[mds$Season=="Wi" & mds$Nacid=="RNA" & mds$Environment=="Coastal"] = 15
  
  points(nmds.s,display="sites",col=nmdscol, pch=pch,cex=.7,lwd=1)
  #text(nmds.s,display="sites",col=nmdscol, cex=0.2, labels=row.names(mds), pos=1)
  legend("bottomright",box.lwd=0.5, pch=c(1,5,0,2,3),
         legend=c("Estuarine Wi","Estuarine Su","Coastal","Offshore","Port"),
         ncol=1,cex=.6)
  legend("topright",box.lwd=0.5, pch=c(1),legend=sort(unique(mds$Estuary)),col=colU,
         ncol=1,cex=.6)
  
  ordihull(nmds.s,groups=mds$Environment,show.groups=c("Coastal","Offshore"),
           col="darkgrey")
  

dev.off()

for (i in c(1:12)) {
  print(efs[[i]])
}

# ----- NMDS with only est. -----------
colE = rb[as.numeric(mdse$Estuary)] 

nmds.e = metaMDS(otus.s.ra.f[mds$Environment=="Estuarine",])
# Run 20 stress 0.1525545 


efse = list(envfit(nmds.e,md1e),envfit(nmds.e,md2e, na.rm=T),
           envfit(nmds.e,md3e, na.rm=T),envfit(nmds.e,md6e, na.rm=T),
           envfit(nmds.e,md7e, na.rm=T),envfit(nmds.e,md8e, na.rm=T),
           envfit(nmds.e,md10e, na.rm=T),envfit(nmds.e,md12e, na.rm=T))

svg("img/sed/NMDS_w_envfit_est_only.svg")

ordiplot(nmds.e,type="none")

for (i in c(1:8)){
  plot(efse[[i]],col="purple",p.max=.001,cex=.5)
}
pch=as.numeric(mdse$Environment)-1
pch[mdse$Season=="Su" & mdse$Nacid=="DNA"] = 5
pch[mdse$Season=="Wi" & mdse$Nacid=="RNA"] = 16
pch[mdse$Season=="Su" & mdse$Nacid=="RNA"] = 18

points(nmds.e,display="sites",col=colE, pch=pch,cex=.7,lwd=1)
text(nmds.e,display="sites",col=colE, cex=0.2, labels=mdse$Station, pos=1)

legend("bottomright",box.lwd=0.5, pch=c(1,5,16,18),
       legend=c("Winter DNA","Summer DNA","Winter RNA","Summer RNA"),
       ncol=1,cex=.6)
       
legend("topright",box.lwd=0.5, pch=c(1),legend=sort(unique(mds$Estuary)),col=colU,
       ncol=1,cex=.6)

dev.off()

for (i in c(1:8)) {
  print(efse[[i]])
}

# ------ Adonis model --------------------

# After AMBI and mgAMBI, redox and PI are the strongest correlated.

otus.m = otus.s.ra.f[!is.na(mds$Redox) & !is.na(mds$OM) & !is.na(mds$PI)& !is.na(mds$Salinity),]
mds.m =  mds[!is.na(mds$Redox) & !is.na(mds$OM) & !is.na(mds$PI) & !is.na(mds$Salinity),]
dim(otus.m) # 287 samples included
adonis(otus.m~Salinity+Redox+OM+PI,data=mds.m)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Salinity    1     7.011  7.0108  28.175 0.07939  0.001 ***
#   Redox       1     8.504  8.5040  34.175 0.09629  0.001 ***
#   OM          1     1.184  1.1841   4.758 0.01341  0.001 ***
#   PI          1     1.443  1.4432   5.800 0.01634  0.001 ***
#   Residuals 282    70.171  0.2488         0.79457           
# Total     286    88.313                 1.00000              



# ------ TAXONOMIC BARCHARTS --------

# Read taxonomic data
taxa.as = read.delim(paste(otudir, "All_Assignments.tsv.gz", sep="/"), header=T, row.names=3)
taxa.ass = as.data.frame(t(taxa.as[,-c(1:2)]))

# Verify correspondence to sample names
table(row.names(mds)==row.names(taxa.ass)) # yes

# Calculate relateive abundances and remove rare taxa
taxa.ass.ra = decostand(taxa.ass,method="total")
taxa.ass.ra.f = dropRareByMaxAbundance(taxa.ass.ra, 10/min(rowSums(taxa.ass)))
taxa.ass.ra.f = decostand(taxa.ass.ra.f, method="total")
dim(taxa.ass.ra.f) # 988 taxa out of 3313

grouping_info<-data.frame(row.names=row.names(mds), mds$Estuary)

# Assignments graph
ta = taxa.ass.ra[mds$Nacid=="DNA",]
table(row.names(mds.dna) == row.names(ta))
row.names(ta)=gsub("_PS","",row.names(ta))
row.names(ta)=gsub("_S_","_",row.names(ta))
grouping_info<-data.frame(row.names=row.names(ta), mds.dna$Estuary)
pdf(paste("img/sed/Taxon_assignments_Sed_DNA.pdf",sep=""),height=7,width=35)
taxaplot(30,grouping_info,ta)
dev.off()



