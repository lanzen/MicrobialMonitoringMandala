# ----- READ DATA -----

setwd("~/projects/uGAMBI_RC/MBC_2019/")
otudir = "SWARM_20190503/CREST_LULU_Final"
require(vegan)
require(vioplot)
source('R/filtering.R')
source('R/diversity.r')
source('R/correlationTests.r')
source('R/taxaplot.R')

mdi = read.table(file="metadata/metadata_final.tsv",header=T,row.names=1,sep="\t")
mdi$intertidal = mdi$PointDepthBM==0
mdi$Date = as.Date(mdi$Date)

otu.file = paste(otudir, "SWARM_table_final.tsv.gz", sep="/")
otus.t = read.delim(otu.file,row.names=1,header=T,sep="\t")
dim(otus.t) #182,502 OTUs, 425 samples

otus.t = otus.t[,order(names(otus.t))]
taxonomy=otus.t$classification
otus = as.data.frame(t(otus.t[,-1]))

mds = droplevels(mdi[mdi$Type=="S",])

mds$AMBI2 = (mds$AMBI+mds$microgAMBI)/2

otus.s = otus[mdi$Type=="S",]

otus.s = otus.s[,colSums(otus.s)>0]
taxonomy.s = taxonomy[colSums(otus.s)>0]

divTemp = read.csv("Diversity_stats/Diversity_classified.csv", row.names=1, header=T)
div.s = divTemp[mdi$Type=="S",c(1,5:7,9:10)]

# Are all samples present?
table(row.names(mds) %in% rownames(otus.s))
# and do they match exactly?
table(row.names(otus.s) == row.names(mds))
table(row.names(div.s) == row.names(mds))

dim(otus.s) # 164,375 OTUs, 251 samples
sum(otus.s) # 8,843,043 reads

summary(div.s)
# 10-138t reads, Rarefied richness = 1144 -- 4623, H' = 5 -- 8.3, J' = 0.63 -- 0.91

## ---- Calculate "partial PIs" for metals and hydrocarbons ----

mds$piMetals = 0
mds$piHC = 0
mds$piOM = NA
mds$piOM5yAvg = NA
metalLimits = data.frame(row.names=c("Zn","Pb","Hg","Cd", "Cr", "Cu", "Ni"),
                         limits = c(249,78,.53,1,39,55,23))
hcLimits = data.frame(row.names=c("TotalPAH","TotalPCB"), 
                      limits = c(1607,24.6))
omLimit = 2.0

geo = c(0,.5,1,2,4,8)
for (s in row.names(mds)){
  metalData = mds[s,row.names(metalLimits)]
  hcData = mds[s,row.names(hcLimits)]
  hasMData = !is.na(metalData)
  hasHCData = !is.na(hcData)
  km = 1/sum(hasMData)
  khc = 1/sum(hasHCData)
  if (sum(hasMData)==0) mds[s,"piMetals"] = NA
  if (sum(hasHCData)==0) mds[s,"piHC"] = NA
  
  for (i in c(1:length(metalData))[hasMData]){
    v = metalData[,i]
    l = metalLimits[i,]
    for (j in c(1:5)){
      if (v<l*geo[j+1]){
        mds[s,"piMetals"] = mds[s,"piMetals"] + km*(v-l*geo[j])/(l*geo[j+1] - l*geo[j])
        break
      } else {
        mds[s,"piMetals"] = mds[s,"piMetals"] + km
      }
    }
  }
  
  for (i in c(1:length(hcData))[hasHCData]){
    v = hcData[,i]
    l = hcLimits[i,]
    for (j in c(1:5)){
      if (v<l*geo[j+1]){
        mds[s,"piHC"] = mds[s,"piHC"] + khc*(v-l*geo[j])/(l*geo[j+1] - l*geo[j])
        break
      } else {
        mds[s,"piHC"] = mds[s,"piHC"] + khc
      }
    }
  }
  if (!is.na(mds[s,"OM"])){
    om = mds[s,"OM"]
    j=1
    piOM = 0
    while(om>omLimit*geo[j] & j<6){
      piOM = piOM + min(1,(om-omLimit*geo[j])/(omLimit*geo[j+1] - omLimit*geo[j]))  
      j=j+1
    }
    mds[s,"piOM"] = piOM
  }
  if (!is.na(mds[s,"OM5yAvg"])){
    om = mds[s,"OM5yAvg"]
    j=1
    piOM5 = 0
    while(om>omLimit*geo[j] & j<6){
      piOM5 = piOM5 + min(1,(om-omLimit*geo[j])/(omLimit*geo[j+1] - omLimit*geo[j]))  
      j=j+1
    }
    mds[s,"piOM5yAvg"] = piOM
  }
}

## ---- Metadata subsets ------

summary(mds)
# 6 from run 2018-09, 113 from 2019, 108 from 2018-12 and 24 from Aylagas 2016 / ItsasTk
# 14 offshore, 14 coastal, 219 estuarine, 4 port
# 39 Summer, 212 winter
# 177 with AMBI and 165 with MAMBI (not ItasasTK) (Summer, 2019, port and SVB missing)
# 227 w AMBI last 5y average, all but port and SVB
# 222 DNA + 29 RNA, Duplicates from 5 samples (technical), one with RNA only
# 48 from Nervion (most samples estuary)
# 146 intertidal and 97 that are always covered

summary(as.factor(mds$Year))
# 2013 2016 2017 2018 2019 
# 61   30   76   49   35 

# Parameters of factor type
mdf = mds[,c("Nacid","Run","Primer","Environment","Estuary","Station","Season",
             "Year","SalClass", "BC.x", "BC.5y","intertidal")] #factors

# All continous parameters
mdc = mds[,-c(1:11,53:54,62:75)]

# Continous parameters available for *all* samples (for NMDS and similar)
md1 = mds[,c("Year","microgAMBI")]

# Parameters grp 2. (all but port and heavily contaminated, n=8 missing)
md2 = mds[,c("AMBI.x","BC.x","Salinity","SalClass","OM5yAvg","OM10yAvg",
              "MaxDepthPM","PointDepthBM")]

# Parameters grp 3 (n=74 missing)
md3 = mds[,c("AMBI","AMBI2")]

# Parameters grp 4 (n=25 missing)
md4 = mds[,c("AMBI.5y","MAMBI.5y")]

# Parameters grp 5 (n=86 missing)
md5= mds[,c("MAMBI","BC.5y")]

# Parameters grp 6 (n=51 missing)
md6= mds[,c("Redox","Cu","Pb","Ni","Cr","Zn","Hg","Cd","piMetals")]

# Parameters grp 7 (n=67 missing)
md7= mds[,c("Mn","Fe")]

# Parameters grp 8 (n=35 missing)
md8= mds[,c("OM","Gravel","Sand","Mud","piOM")]

# Parameters grp 9 (n=21 missing)
md9= mds[,c("PI5yAvg","PI10yAvg","Redox5yAvg","Redox10yAvg")]

# Parameters grp 10 (n=47 missing)
md10= mds[,c("PI","TotalPCB","TotalPAH","TotalDDT","TotalHCH","piHC")]

# Parameters grp 11 (n=22 missing)
md11= mds[,c("WaterO2","WaterO2")]

# Parameters grp 12 (n=41 missing)
md12= mds[,c("Q1SiteSal","Q2SiteSal","Q3SiteSal")]

# Calculate relative abundance and Hellinger transformed and filter rare
otus.s.ra.all = decostand(otus.s, method="total")
otus.s.ra.f1 = dropRareByMaxAbundance(otus.s.ra.all, 10/min(rowSums(otus.s)))
otus.s.ra.f = decostand(otus.s.ra.f1, method="total")
dim(otus.s.ra.f) # 4741 otus of >180k retained

taxonomy.s.ra.f = as.character(otus.t[names(otus.s.ra.f),"classification"])

otus.s.H = decostand(otus.s.ra.f1, method="hell")
otus.s.ra.bc = vegdist(otus.s.ra.f)

## Choose variables to train on for ML and TITAN/Splines:
trainOn = c("AMBI", "MAMBI", "microgAMBI","PI", "PI5yAvg", "piMetals", "piHC", "piOM", "piOM5yAvg", "Redox","Redox5yAvg")

# Graph colours
rb = rainbow(17)
colour = rb[as.numeric(mds$Estuary)] 
colU = rb[sort(as.numeric(unique(mds$Estuary)))]
# make None = black
colour[mds$Estuary == "None"] = "black"
colU[sort(unique(mds$Estuary)) == "None"] = "black"

# -------- Check distribution of rare and dominant OTUs ------------------

otus.lessAb = otus.s[,!(names(otus.s) %in% names(otus.s.ra.f1))]

otus.rarest = otus.s[,colMeans(otus.s.ra.all)<median(colMeans(otus.s.ra.all))]

AbCrossSamplDist = colSums(decostand(otus.s.ra.f,method="pa"))
AllCrossSamplDist = colSums(decostand(otus.s,method="pa"))
LessABCrossSDist = colSums(decostand(otus.lessAb,method="pa"))
rarestCrossSDist = colSums(decostand(otus.rarest,method="pa"))

summary(AbCrossSamplDist)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   28.00   61.00   75.27  114.00  250.00

summary(AllCrossSamplDist)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   2.000   3.000   9.833   8.000 250.000 

summary(rarestCrossSDist)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.00   1.00   2.00   1.88   2.000   7.00 

pdf("img/sed/no_of_samples_where_OTU_present_distro.pdf")
vioplot(AbCrossSamplDist,AllCrossSamplDist,LessABCrossSDist,
        names=c("Abundant OTUs","All OTUs","Below median OTUs"))
dev.off()

# hist(AbCrossSamplDist,breaks=100)
# hist(AllCrossSamplDist,breaks=100)

wilcox.test(AbCrossSamplDist,AllCrossSamplDist)
# p<2E-16

# -------- Make subsets for removing replicates, DNA or RNA only ----------

# Remove all B replicates (for parameter correlation, n=247)
mds.nr = mds[mds$Rep!="B",]
div.s.nr = div.s[mds$Rep!="B",]
otus.s.ra.f.nr = otus.s.ra.f[mds$Rep!="B",]

# Remove all B replicates and RNA (for parameter correlation, n=217)
mds.nr.dna = mds.nr[mds.nr$Nacid=="DNA",]
mdf.nr.dna = mds.nr.dna[,c("Nacid","Run","Primer","Environment","Estuary","Station","Season",
             "Year","SalClass", "BC.x", "BC.5y","intertidal")]
div.s.nr.dna = div.s.nr[mds.nr$Nacid=="DNA",]
otus.s.ra.f.nr.dna = otus.s.ra.f.nr[mds.nr$Nacid=="DNA",]

# Technical replicates only subset (n=5)
mds.repsA = mds[mds$Rep=="A",]
mds.repsB = mds[mds$Rep=="B",]
div.s.repsA = div.s[mds$Rep=="A",]
div.s.repsB = div.s[mds$Rep=="B",]
otus.s.ra.f.repsA = otus.s.ra.f[mds$Rep=="A",]
otus.s.ra.f.repsB = otus.s.ra.f[mds$Rep=="B",]
otus.s.ra.f.repsBoth = rbind(otus.s.ra.f.repsA,otus.s.ra.f.repsB)
bcRepDiff = diag(as.matrix(vegdist(otus.s.ra.f.repsBoth))[c(1:5),c(6:10)])
summary(bcRepDiff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1064  0.1529  0.1551  0.1513  0.1682  0.1737 

# DNA vs RNA only subset (for samples were both were sequenced, n=28)
mds.dnaComp = mds[mds$DNA_RNA=="D",]
mds.rnaComp = mds[mds$DNA_RNA=="R",]
div.s.dnaComp = div.s[mds$DNA_RNA=="D",]
div.s.rnaComp = div.s[mds$DNA_RNA=="R",]
otus.s.ra.f.dnaComp = otus.s.ra.f[mds$DNA_RNA=="D",]
otus.s.ra.f.rnaComp = otus.s.ra.f[mds$DNA_RNA=="R",]
otus.s.ra.f.RNAvDNA = rbind(otus.s.ra.f.dnaComp,otus.s.ra.f.rnaComp)
bcRNAvDNA_Diff = diag(as.matrix(vegdist(otus.s.ra.f.RNAvDNA))[c(1:28),c(29:56)])
summary(bcRNAvDNA_Diff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2658  0.3358  0.3886  0.3846  0.4153  0.5457 

interSampleBC = as.numeric(unlist(vegdist(otus.s.ra.f.nr.dna)))
summary(interSampleBC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1807  0.6738  0.7947  0.7705  0.8876  0.9948

wilcox.test(bcRepDiff,interSampleBC,alternative="less") # p=5E-5, replicates are more similar
wilcox.test(bcRNAvDNA_Diff,interSampleBC,alternative="less") # p<2E-16, DNA-RNA are more similar
wilcox.test(bcRNAvDNA_Diff,bcRepDiff) # p=8.4E-6, DNA-RNA are more different than tech reps

pdf("img/sed/BCDifs_viol_replication_DNA_RNA.pdf")
vioplot(bcRepDiff,bcRNAvDNA_Diff,interSampleBC,names=c("Reps","RNA v DNA","Inter-sample"))
dev.off()

# -------- MicrogAMBI vs. AMBI, Variation between replicates, DNA vs RNA -----

repCol = colour[mds$Rep=="A"]
lmMgA = lm(microgAMBI~AMBI,data=mds)
summary(lmMgA) # p<2E-16, adj R2=.44, 

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
       pch=c(rep(1,252),2,1,3),cex=.6,ncol=1)

# text(mds$microgAMBI~mds$AMBI, col=colour, pch=as.numeric(mds$Environment),cex=.3,
#      label=mds$Station, pos=1)
# text(mds$microgAMBI~mds$AMBI, col=colour, pch=as.numeric(mds$Environment),cex=.3,
#      label=mds$Year, pos=3)

dev.off()

techRepDiff= abs(mds.repsA$microgAMBI - mds.repsB$microgAMBI)
summary(techRepDiff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.006766 0.012220 0.025280 0.026070 0.034040 0.052020 

lmMgD = lm(microgAMBI~AMBI,data=mds.nr.dna)
summary(lmMgD) #Adj R2=.44, mg = .71 + AMBI*.50 (almost identical to including all)
lmMgR = lm(microgAMBI~AMBI,data=mds.rnaComp)
summary(lmMgR) #Adj R2=-.02, p=.47, mg=2.5+AMBI*.12 (very poor)
lmMgDComp = lm(microgAMBI~AMBI,data=mds.dnaComp)
summary(lmMgDComp) #Adj R2=0, p=.34, mg=2.3+0.13*AMBI (v poor)


pdf("img/sed/microgAMBI_vs_AMBI_RNA_vs_DNA_all.pdf",width=6,height = 6)
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
#abline(lmMgD,col="red",lty=2) # too similar
legend("topleft",pch=c(1,16),col=c(2,4),legend=c("DNA","RNA"))
dev.off()

# text(mds$microgAMBI~mds$AMBI, col=colour, pch=as.numeric(mds$Environment),cex=.3,
#      label=mds$Station, pos=1)
# text(mds$microgAMBI~mds$AMBI, col=colour, pch=as.numeric(mds$Environment),cex=.3,
#      label=mds$Year, pos=3)

dev.off()

# Is there a difference between how mgAMBI is predicted btw DNA v RNA?
RNAmgDiff= mds.rnaComp$microgAMBI - mds.dnaComp$microgAMBI
summary(RNAmgDiff)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.77110 -0.03455  0.23320  0.14560  0.35690  0.81530 
wilcox.test(RNAmgDiff) #p=0.01
# --> Yes, in average RNA predicts higher mgAMBI than DNA!

# Does RNA or DNA predicted mgAMBI better resemble AMBI?
RNA_DNA_mgAMBI_v_AMBI = data.frame(RNADiff = mds.rnaComp$microgAMBI - mds.rnaComp$AMBI,
                                   DNADiff = mds.dnaComp$microgAMBI - mds.dnaComp$AMBI,
                                   RNADiffAbs = abs(mds.rnaComp$microgAMBI - mds.rnaComp$AMBI),
                                   DNADiffAbs = abs(mds.dnaComp$microgAMBI - mds.dnaComp$AMBI))

summary(RNA_DNA_mgAMBI_v_AMBI)
# RNADiff           DNADiff           RNADiffAbs      DNADiffAbs     
# Min.   :-3.2146   Min.   :-3.67104   Min.   :0.162   Min.   :0.01852  
# 1st Qu.:-1.7074   1st Qu.:-1.65294   1st Qu.:0.740   1st Qu.:0.68198  
# Median :-0.9895   Median :-0.93964   Median :1.033   Median :1.04540  
# Mean   :-0.8394   Mean   :-1.03163   Mean   :1.281   Mean   :1.25361  
# 3rd Qu.: 0.2783   3rd Qu.:-0.02087   3rd Qu.:1.707   3rd Qu.:1.65294  
# Max.   : 1.5424   Max.   : 1.20817   Max.   :3.215   Max.   :3.67104  
# NA's   :3         NA's   :3          NA's   :3       NA's   :3        

wilcox.test(RNA_DNA_mgAMBI_v_AMBI$RNADiff,
            RNA_DNA_mgAMBI_v_AMBI$DNADiff,alternative="less")
wilcox.test(RNA_DNA_mgAMBI_v_AMBI$RNADiffAbs,
            RNA_DNA_mgAMBI_v_AMBI$DNADiffAbs,alternative="less") 

# --> Yes, both underpredict but RNA less, so closer to AMBI, but far from significant (p=.6)


# -------- Parameter correlations  ---------  

# Sanity check of run and primer, which makes little sense because 
# samples remaining from Eva's run are mainly coastal with low impact
# printANOVA(mdf[,c(2:3)],mds[,-c(1:11,53:54,62:74)],a=.05/48)

printANOVA(mdf.nr.dna[,-c(1:3,5:6)],mds.nr.dna[,-c(1:11,53:54,62:74)],a=.05/(47*7))
# Coastal and offshore have lower microgAMBI than estuaries (p=3E-7)
# Same for AMBI (p=2E-14), its 5y average (p<.02) and AMBI.x (p=4E-16)
# Estuaries have lower redox v coast (p=3E-6) and more mud (p=3E-5)
# Redox is lower in summer (below 0, p=6E-5) and water O2 too ("0")
# Salinity class follows clear trend of lower AMBI
# Intertidal have lower metal conc. esp Hg and lower PCB and OM, lower PI and higher redox
# This may be due to the majorirty of estuaries having intertidal sampling

printAllvsAll(mds.nr.dna[,-c(1:11,53:54,62:75)],.05/1128)
# Year on year trends unreliable due to sampling bias
# AMBI ~ mgAMBI correlates better than mgAMBI ~ AMBI.5y (R2=0.49, 0.45)
# MAMBI ~ mgAMBI  much weaker (R2=.10, p=8E-5 or 0.14, 5E-8 for MAMBI.5y)
# AMBI ~ MAMBI: R2=.48
# AMBI.5y ~ MAMBI.5y; R2=.56

# mgAMBI ~ Redox: p=4E-12, R2=.25 (better than AMBI)
# AMBI ~ Redox: p=5E-9, R2=.22
# mgAMBI ~ Redox.5y: p<2E-16, R2=.50 (much better than yearly, more integrative than AMBI)
# AMBI ~ Redox.5y: p=2E-13, R2=.33 (also better than yearly though)
# mgAMBI ~ Redox.10y: p<2E-16, R2=.49
# AMBI ~ Redox.10y: p=3E-14, R2=.35
# AMBI.x ~ Redox: R2=.41
# AMBI.x ~ Redox.5y: R2=.58
# AMBI.x ~ Redox.10y: R2=.60 (best)
# AMBI2 ~ Redox: R2=.29 (better than either AMBI or mgAMBI)
# AMBI2 ~ Redox.5y: R2=.47
# AMBI2 ~ Redox.10y: R2=.50

# mgAMBI ~ OM: p=6E-7, R2=.13 (comparable to AMBI) (PI_OM: 2E-11, R2=.22)
# AMBI ~ OM: p=4E-6, R2=.13 (PI_OM: 7E-7, R2=.15)
# mgAMBI ~ OM.5y: p=9E-10, R2=.16 (comparable to AMBI, better than yearly)
# AMBI ~ OM.5y: p=8E-8, R2=.17
# mgAMBI ~ OM.10y: p=1E-11, R2=.20 (better than yearly / mgAMBI integrative, like AMBI)
# AMBI ~ OM.10y: p=7E-10, R2=.23
# AMBI.5y ~ OM.5y: R2=.46 (indicates very unreliable OM but with strong influence)
# AMBI.x ~ OM: R2=.212
# AMBI.x ~ OM.5y: R2=.31 
# AMBI.x ~ OM.10y: R2=.36 (best)
# MAMBI ~ OM.10y: o=9E-8, R2=.16 (not as good as AMBI)
# AMBI2 ~ OM: p=4E-7, R2=.13 (better)
# AMBI2 ~ OM.5y: p=3E-9, R2=.18 (better)
# AMBI2 ~ OM.10y: p=8E-12, R2=.23 (better)

# mgABMI ~ mud: p=2E-10, R2=0.20 (comparable to AMBI)
# AMBI ~ Mud: p=4E-8, R2=.18
# AMBI.x ~ mud: R2=.25
# AMBI2 ~ mud: p=5E-12, R2=.24 (best)

# mgAMBI ~ PI: p=9E-10, R2=0.19 (better than AMBI but 5y AMBI average best)
# AMBI ~ PI: p=3E-6, R2=.14
# mgAMBI ~ PI.5y: p=2E-16, R2=.28 (better than yearly / integrative, much better than AMBI)
# AMBI ~ PI.5y: p=9E-7, R2=.16 (not better than yearly like mgAMBI)
# mgAMBI ~ PI.10y: p=3E-15, R2=.27
# AMBI ~ PI.10y: p=8E-8, R2=.19
# AMBI.5y ~ PI.5y: p=9E-12, R=.21
# AMBI.x ~ PI: R2=.33
# AMBI.x ~ PI.5y: R2=.37
# AMBI.x ~ PI.10y: R2=.39 (best)
# AMBI2 ~ PI: p=5E-10, R2=.21 (better than mgAMBI or AMBI)
# AMBI2 ~ PI.5y: R2=.23 (not as impressive), 10y: R2=.26 (slightly better)

# mgAMBI ~ PI_Metals: p=2E-6, R2=.12
# AMBI5y ~ PI_Metals: 1.5E-8, R2=.18

# mgAMBI ~ PI_HC: p=4E-6m R2=.11

# mgAMBI ~ salinity (neg): R2=.13, p=7E-8 (even stronger for AMBI, R2=.19)

# Very interesting that mgAMBI appears to correlate stronger with 5y averages for PI, OM 
# and redox. Though AMBI does also, except perhaps weaker and not for PI.

# Redox ~ OM (neg): p<2E-16, R2=.37 (less for redox 5 or 10y) (better fit with PI_OM:R2=.48)
# Redox ~ Mud (neg): p<2E-16, R2=.45 (5y Redox .38)
# Redox ~ PI (neg):  p<2E-16, R2=.40 (5y and 10y comparable) 

# Mud ~ OM: p=7E-13, R2=.24 (much stronger, R2=.41 with PI_OM)
# PI ~ OM: p=2E-14, R2=.29
# Mud ~ PI: p=4E-13, R2=.26 (5y similar)
# Metals, PCBs etc to each other


# ------- Alphadiv correlation -----

printANOVA(mdf.nr.dna, div.s.nr.dna, .05)

# Higher Chao1 with new primers
# RNA: BC.x == 4 has lower rarefied S and H' and J' vs 2 (highest, p=0)
# Same for BC.5y but with BC==5 too (and 4 for H' and J')

printVS(div.s.nr.dna, mds.nr.dna[,-c(1:11,53:54,62:75)], a=.05/49)
# Trends of more reads in higher salinity, lower impact
#RS ~ Pb, Cd
#             p       Adj R2 (Adj R2 with RNA)
# H'~mgAMBI   2E-8   .13    (0.17)
# (H'~AMBI5y   9E-5    .06)
# H'~AMBIx    2E-6    .09   (.12) (but local maximum around 2?)
# H'~Cd       1E-3    .0.55 (.51)
# H~PIMetals  7E-5    .085
# J'~mgAMBI   4E-12   .20   (.23)
# J'~AMBI5y   2E-5    .087  (.11)
# J'~AMBIx    6E-10   .16   (.18)
# J'~Cu       9E-4    .059  (.059)
# J~PI        8E-5    .082  (.090)
# J~PIMetals  6E-5    .088
# J~PIHC      2E-4    .071
# J~PI5y      9E-5    .070  (.079)
# J~PI10y     9E-4    .049  (.059)
# J~Redox5y   3E-5    .082  (.090)
# J~Redox10y  1E-4    .069  (.089)
# J~OM10y     3E-4    .058  (.058)
# D~mgAMBI    3E-9    .15   (.15)
# D~OM10y     5E-4    .045
# Chao~AMBI   1E-5    .10 (positive!)

repCol = colour[mds$Rep=="A"]

m1 = glm(div.s$J~mds$microgAMBI)

pdf("img/sed/microgAMBI_vs_J.pdf",width=6,height = 6)
plot(div.s$J ~ mds$microgAMBI, col=colour, pch=as.numeric(mds$Environment),cex=.6,
     xlab="microgAMBI",ylab="Pielou evenness")
for (i in c(1:dim(mds.repsA)[1])){
  lines(x=c(mds.repsA$microgAMBI[i],mds.repsB$microgAMBI[i]),
        y=c(div.s.repsA$J[i],div.s.repsB$J[i]),
        col=repCol[i])
}
abline(lmMgJ,col="grey",lty=2)

legend("bottomleft",col=c(colU,"white",rep("black",3)), 
       legend=(c(as.vector(sort(unique(mdi$Estuary))),"","Estuary","Coast","Offshore")),
       pch=c(rep(1,252),2,1,3),cex=.6,ncol=3)

dev.off()

m1 = glm(div.s$J~mds$PI)

pdf("img/sed/PI_vs_evenness.pdf",width=6,height = 6)
plot(div.s$J ~ mds$PI, col=colour, pch=as.numeric(mds$Environment),cex=.6,
     xlab="PI",ylab="Pielou evenness")
for (i in c(1:dim(mds.repsA)[1])){
  lines(x=c(mds.repsA$PI[i],mds.repsB$PI[i]),
        y=c(div.s.repsA$J[i],div.s.repsB$J[i]),
        col=repCol[i])
}
abline(lmMgJ,col="grey",lty=2)
dev.off()

# ------- Alphadiv correlation, DNA v RNA -----

# Are ther any differences in correlations w regard to DNA v RNA diversity and paramters

printANOVA(mds.dnaComp[,c("Environment","Estuary","Season",
                           "BC.x", "BC.5y")], div.s.dnaComp[,-1], .05)

printANOVA(mds.rnaComp[,c("Environment","Estuary","Season",
                           "BC.x", "BC.5y")], div.s.rnaComp[,-1], .05)

# RNA: BC.x == 2 has higher rarefied S and H' and J' (p<.007) and 1-D (p=.04)

printVS(div.s.dnaComp[,-1], mds.dnaComp[,-c(1:11,53:54,62:73)], a=.01)
printVS(div.s.rnaComp[,-1], mds.rnaComp[,-c(1:11,53:54,62:73)], a=.01)

#              DNA p   Adj.R2 |  RNA p   Adj R2
# RS~mgAMBI                      2E-4    .50
# RS~AMBIx                       .004    .24
# RS~Ni                          .002    .29
# RS~Cr                          .01     .18
# J'~mgAMBI    .04     .12       2E-5    .49
# J'~AMBI5y                      .02     .17
# J'~AMBIx                       .002    .28
# J'~Ni                          5E-4    .35
# J'~Cr                          .01     .18
# J'~PI                          .04     .11
# J'~PI5y      .004    .25
# J'~PI10y                       .03     .14 
# J'~Redox5y   .006    .23       .005    .24
# J'~Redox10y  .02     .16       .002    .27
# (1-D)~mgAMBI .03     .13       7E-4    .34
# (1-D)~AMBIx  .02     .15       .02     .15
# (1-D)~Ni     .03     .14       .005    .24
# (1-D)~Zn     .01     .20
# (1-D)~Cd     .03     .14
# (1-D)~PCB    .04     .12
# (1-D)~PI     .02     .16
# (1-D)~PI5y   .07     .08

# (H'~mgAMBI                      2E-6    .56
# H~'AMBIx                       .01     .20
# ...)

div.s.RNAvDNA = (div.s.rnaComp/div.s.dnaComp)[,-1]
summary(div.s.RNAvDNA)
# Rarefied.richness       H                J             Simpson           Chao1       
# Min.   :0.7336    Min.   :0.8757   Min.   :0.9278   Min.   :0.9933   Min.   :0.4688  
# 1st Qu.:0.8632    1st Qu.:0.9443   1st Qu.:0.9760   1st Qu.:0.9981   1st Qu.:0.6541  
# Median :0.9169    Median :0.9814   Median :0.9918   Median :0.9996   Median :0.7247  
# Mean   :0.9197    Mean   :0.9725   Mean   :0.9948   Mean   :0.9994   Mean   :0.7230  
# 3rd Qu.:0.9568    3rd Qu.:0.9950   3rd Qu.:1.0188   3rd Qu.:1.0004   3rd Qu.:0.8048  
# Max.   :1.1341    Max.   :1.0296   Max.   :1.0595   Max.   :1.0064   Max.   :1.0578  

div.s.DNAvRNADiff = (div.s.dnaComp-div.s.rnaComp)[,-1]
summary(div.s.DNAvRNADiff)
# Rarefied.richness       H                 J                Simpson          
# Min.   :-485.3    Min.   :-0.2101   Min.   :-0.051078   Min.   :-0.0063278  
# 1st Qu.: 138.7    1st Qu.: 0.0355   1st Qu.:-0.015059   1st Qu.:-0.0004351  
# Median : 265.6    Median : 0.1336   Median : 0.007043   Median : 0.0004382  
# Mean   : 274.3    Mean   : 0.2016   Mean   : 0.004510   Mean   : 0.0005616  
# 3rd Qu.: 494.9    3rd Qu.: 0.4122   3rd Qu.: 0.020778   3rd Qu.: 0.0019047  
# Max.   : 942.8    Max.   : 0.9300   Max.   : 0.058984   Max.   : 0.0065838  
# Chao1        
# Min.   : -592.1  
# 1st Qu.: 2595.8  
# Median : 3816.3  
# Mean   : 4751.7  
# 3rd Qu.: 5533.1  
# Max.   :14125.3

lmDRRS=lm(div.s.dnaComp$Rarefied.richness~div.s.rnaComp$Rarefied.richness)
summary(lmDRRS) #p=3E-6, Adj R2.55
º

lmRSQMud = lm(div.s.RNAvDNA$Rarefied.richness~mds.dnaComp$Mud)
summary(lmRSQMud)
plot(div.s.RNAvDNA$Rarefied.richness~mds.dnaComp$Mud)

lmRSQR10y = lm(div.s.RNAvDNA$Rarefied.richness~mds.dnaComp$Redox10yAvg)
summary(lmRSQR10y) # p=.001, Adj R2=.31
pdf("img/sed/RNAvDNA_RS_vs_Redox10y.pdf",width=4.8,height=5)
plot(div.s.RNAvDNA$Rarefied.richness~mds.dnaComp$Redox10yAvg, 
     xlab="Redox (10y arithmetic mean)",ylab="RNA/DNA rarefied richness")
abline(lmRSQR10y,col="grey",lty=2)
dev.off()


pdf("img/sed/microgAMBI_vs_J_RNA_vs_DNA_all.pdf",width=6,height = 6)
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


pdf("img/sed/PI_vs_J_RNA_vs_DNA_all.pdf",width=6,height = 6)
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

# ------- ANOSIM and adonis  --------

anosim(otus.s.ra.f.nr.dna, grouping=mds.nr.dna$Station)
# R=.72 ***

anosim(otus.s.ra.f.nr.dna, grouping=mds.nr.dna$Estuary)
# R=.28***

anosim(otus.s.ra.f.nr.dna, grouping=mds.nr.dna$Environment)
# R=.39 ***

anosim(otus.s.ra.f.nr.dna, grouping=mds.nr.dna$Season)
# NS

anosim(vegdist(otus.s.ra.f.nr.dna[!is.na(mds.nr.dna$SalClass),]), 
       grouping=mds.nr.dna$SalClass[!is.na(mds.nr.dna$SalClass)])
#R=.23***

adonis(otus.s.ra.f.nr.dna~Environment+Season+Station+Year,data=mds.nr.dna)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Environment   3     8.357 2.78569 19.9113 0.12586  0.001 ***
#   Season        1     1.765 1.76497 12.6155 0.02658  0.001 ***
#   Station      63    35.143 0.55782  3.9872 0.52928  0.001 ***
#   Year          1     0.426 0.42643  3.0480 0.00642  0.004 ** 
#   Residuals   148    20.706 0.13990         0.31185           
# Total       216    66.397                 1.00000    

# ------- adonis RNA v DNA ---------

# Is there a difference in seasonal difference DNA v RNA?
anosim(otus.s.ra.f.rnaComp,mds.rnaComp$Season)
# R=.53, p=.015
anosim(otus.s.ra.f.dnaComp,mds.dnaComp$Season)
# R=.54, p=.011
# Indicates weaker seasonaility of RNA (contra-intuitive but weak)

# Is there a difference between how DNA v RNA is structured with AMBI?
adonis(otus.s.ra.f.rnaComp[!is.na(mds.rnaComp$AMBI),]~AMBI,data=mds.rnaComp[!is.na(mds.rnaComp$AMBI),])
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# AMBI       1    0.3831 0.38313   1.601 0.06508  0.072 .
# Residuals 23    5.5040 0.23930         0.93492         
# Total     24    5.8871                 1.00000     

adonis(otus.s.ra.f.dnaComp[!is.na(mds.dnaComp$AMBI),]~AMBI,data=mds.dnaComp[!is.na(mds.dnaComp$AMBI),])
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# AMBI       1    0.3231 0.32314    1.38 0.05661  0.145
# Residuals 23    5.3855 0.23415         0.94339       
# Total     24    5.7087                 1.00000

# Indicates better structuring with AMBI for RNA

# Is there a difference between how DNA v RNA is structured with PI?
adonis(otus.s.ra.f.rnaComp[!is.na(mds.rnaComp$PI),]~PI,data=mds.rnaComp[!is.na(mds.rnaComp$PI),])
  #R2=0.101
adonis(otus.s.ra.f.dnaComp[!is.na(mds.dnaComp$PI),]~PI,data=mds.dnaComp[!is.na(mds.dnaComp$PI),])
    #R2=.104

adonis(otus.s.ra.f.rnaComp[!is.na(mds.rnaComp$Redox10yAvg),]~Redox10yAvg,data=mds.rnaComp[!is.na(mds.rnaComp$Redox10yAvg),])
#R2=.116**

adonis(otus.s.ra.f.dnaComp[!is.na(mds.dnaComp$Redox10yAvg),]~Redox10yAvg,data=mds.dnaComp[!is.na(mds.dnaComp$Redox10yAvg),])
#R=.108***

# ----- NMDS -----

nmds.s = metaMDS(otus.s.ra.f)

efs = list(envfit(nmds.s,md1),envfit(nmds.s,md2,na.rm=T),
        envfit(nmds.s,md3,na.rm=T),envfit(nmds.s,md4,na.rm=T),envfit(nmds.s,md5,na.rm=T),
        envfit(nmds.s,md6,na.rm=T),envfit(nmds.s,md7,na.rm=T),envfit(nmds.s,md8,na.rm=T),
        envfit(nmds.s,md9,na.rm=T),envfit(nmds.s,md10,na.rm=T),envfit(nmds.s,md11,na.rm=T),
        envfit(nmds.s,md12,na.rm=T))

svg("img/sed/NMDS_w_envfit.svg")
  
  ordiplot(nmds.s,type="none")
  for (i in c(1:3,6:8, 10:11)){
    plot(efs[[i]],col="purple",p.max=.001,cex=.5)
  }
  nmdscol=colour
  #nmdscol=as.numeric(mds$Nacid)*2
  points(nmds.s,display="sites",col=nmdscol, pch=as.numeric(mds$Environment),cex=.7,lwd=2)
  #text(nmds.s,display="sites",col=nmdscol, cex=0.2, labels=row.names(mds), pos=1)
  legend("bottomright",box.lwd=0.5, pch=c(2,1,3,4),legend=c("Estuarine","Coastal",
                                                            "Offshore","Port"),
         ncol=1,cex=.6)
  legend("topright",box.lwd=0.5, pch=c(1),legend=sort(unique(mds$Estuary)),col=colU,
         ncol=1,cex=.6)
  
  ordihull(nmds.s,groups=mds$Environment,show.groups=c("Coastal","Offshore"),
           col="darkgrey")

dev.off()

for (i in c(1:12)) {
  print(efs[[i]])
}

# ------ Adonis model --------------------

# After AMBI and mgAMBI, redox and PI are the strongest correlated. Even more for 5y averages

otus.m = otus.s.ra.f[!is.na(mds$Redox) & !is.na(mds$OM) & !is.na(mds$PI)& !is.na(mds$Salinity),]
mds.m =  mds[!is.na(mds$Redox) & !is.na(mds$OM) & !is.na(mds$PI) & !is.na(mds$Salinity),]
dim(otus.m) # 196 samples included
adonis(otus.m~Salinity+Redox+OM+PI,data=mds.m)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Salinity    1     4.892  4.8915 19.9977 0.08134  0.001 ***
#   Redox       1     6.555  6.5546 26.7965 0.10899  0.001 ***
#   OM          1     0.975  0.9755  3.9879 0.01622  0.001 ***
#   PI          1     0.999  0.9988  4.0833 0.01661  0.001 ***
#   Residuals 191    46.720  0.2446         0.77685           
# Total     195    60.140                 1.00000             

otus.m = otus.s.ra.f[!is.na(mds$Redox5yAvg) & !is.na(mds$OM5yAvg) & !is.na(mds$PI5yAvg)& !is.na(mds$Salinity),]
mds.m =  mds[!is.na(mds$Redox5yAvg) & !is.na(mds$OM5yAvg) & !is.na(mds$PI5yAvg)& !is.na(mds$Salinity),]
dim(otus.m) # 231 samples included
adonis(otus.m~Salinity+Redox5yAvg+OM5yAvg+PI5yAvg,data=mds.m)

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Salinity     1     5.350  5.3497  22.635 0.07722  0.001 ***
#   Redox5yAvg   1     8.742  8.7415  36.986 0.12618  0.001 ***
#   OM5yAvg      1     0.953  0.9525   4.030 0.01375  0.001 ***
#   PI5yAvg      1     0.821  0.8213   3.475 0.01186  0.001 ***
#   Residuals  226    53.415  0.2363         0.77100           
# Total      230    69.280                 1.00000           

# If only inlcuding those 195 samples with also current values, residuals drop to 76.8%

# ------ TAXONOMIC BARCHARTS --------

# Read taxonomic data
taxa.all = read.delim(paste(otudir, "Relative_Abundance.tsv.gz", sep="/"), header=T, row.names=3)
table(row.names(mdi)==names(taxa.all)[-c(1:2)]) # yes

taxa.as = read.delim(paste(otudir, "All_Assignments.tsv.gz", sep="/"), header=T, row.names=3)

taxa.s = taxa.all[,c(TRUE, TRUE, mdi$Type=="S")]

taxa.ass = as.data.frame(t(taxa.as[,c(FALSE, FALSE, mdi$Type=="S")]))

table(row.names(mds)==names(taxa.s)[-c(1:2)]) # yes
table(row.names(mds)==row.names(taxa.ass)) # yes

taxa.ass.ra = decostand(taxa.ass,method="total")
taxa.ass.ra.f = dropRareByMaxAbundance(taxa.ass.ra, 10/min(rowSums(taxa.ass)))
taxa.ass.ra.f = decostand(taxa.ass.ra.f, method="total")
dim(taxa.ass.ra.f) # 925 taxa out of 3185

grouping_info<-data.frame(row.names=row.names(mds), mds$Estuary)
ranks = data.frame(rank=c("domain","superkingdom","kingdom","phylum","class","order","family","genus","species"),
                   levels=c(3,11,26,26,26,26,26,26,26))

for (i in c(1:dim(ranks)[1])){
  level=toString(ranks$rank[i])
  levelTaxa = as.data.frame(t(taxa.s[taxa.s$Rank==level,-c(1:2)]))
  #row.names(levelTaxa) = row.names(mds)
  print(paste(mean(rowSums(levelTaxa))*100,"% classified at rank",level))
  #line <- readline()      
  pdf(paste("img/sed/16SBarcharts/",level,"_Sed.pdf",sep=""),height=7,width=30)
  taxaplot(ranks$levels[i],grouping_info,levelTaxa)
  dev.off()
}
# [1] "99.9995766846109 % classified at rank domain"
# [1] "99.9995766846109 % classified at rank superkingdom"
# [1] "99.8725983728347 % classified at rank kingdom"
# [1] "99.7075428332733 % classified at rank phylum"
# [1] "96.181818946663 % classified at rank class"
# [1] "83.2638481032804 % classified at rank order"
# [1] "68.5444519647607 % classified at rank family"
# [1] "29.9591219434818 % classified at rank genus"
# [1] "0.187169323624386 % classified at rank species"

# ---
st = as.data.frame(t(taxa.s[taxa.s$Rank=="species",-c(1:2)]))
gt = as.data.frame(t(taxa.s[taxa.s$Rank=="genus",-c(1:2)]))

vs = data.frame(V_maisflavi=st$`Vibrio marisflavi CECT 7928`,
                V_agarivorans=st$`Vibrio agarivorans`,
                Vibrio_Unclass = gt$Vibrio - st$`Vibrio marisflavi CECT 7928` - st$`Vibrio agarivorans`,
                row.names=mds$Sample)
pdf("img/16S_barcharts/Vibrio_sed.pdf",height=8,width=20)
taxaplot(3,grouping_info,vs)
dev.off()

# ------ microgAMBI ~ physchem or time -------

#stations = unique(mds$Station)
summer = mds$Station[mds$Season=="Su" & mds$Estuary!="00 SVB"]
stations = unique(summer)
n = length(stations) #18

sc =as.data.frame(mds)
sc$keep =FALSE
divc = as.data.frame(div.s)
divc$keep = FALSE

sc$colour = colour

sc$Date = as.Date(sc$Date)

stn=as.character(stations[1])
sc[sc$Station==stn,]$keep = TRUE
divc[sc$Station==stn,]$keep = TRUE
st = sc[sc$Station==stn,]
st = st[order(st$Date),]

sc=droplevels(sc[sc$keep,])
divc = divc[divc$keep,]


# ---- Statistics on winter vs summer ----

printANOVA1Factor(sc$Season, divc, a=0.05)
# Chao1 slightly higher in winter, nothing else

printANOVA1Factor(sc$Season, sc[,c(13,15:17,19,21:35,41:52)],a=.05)
# Higher mgAMBI in winter, also higher O2

pdf("img/sed/Seasonal/Wi_v_Su_mgAMBI.pdf")
vioplot(sc$microgAMBI[sc$Season=="Wi"],sc$microgAMBI[sc$Season=="Su"], names=c("Winter","Summer"))
dev.off()

# ---- PI and AMBI ----

stn=as.character(stations[1])
st = sc[sc$Station==stn,]
st = st[order(st$Date),]
plot(st$PI~st$Date,type="l",xlim=range(as.Date("2013-01-01"),as.Date("2018-10-30")),
     ylim=range(0,4.5),col=st$colour, pch=as.numeric(st$Season),xlab="Date",ylab="PI")
for (i in c(2:n)) {
  stn=as.character(stations[i])
  st = sc[sc$Station==stn,]
  if (dim(st)[1]>1){
    st = st[order(st$Date),]
    lines(st$PI~st$Date,type="l",
          col=st$colour,pch=as.numeric(st$Season))
    last = st[dim(st)[1],]
    text(last$Date,last$PI,
         col=last$colour,pos=4,labels = last$Station,cex=.5)
  }
}



