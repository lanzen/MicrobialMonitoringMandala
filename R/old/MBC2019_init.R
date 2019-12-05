setwd("~/projects/uGAMBI_RC/MBC_2019/")
otudir = "SWARM_20190503/CREST_LULU"
require(vegan)

source('R/filtering.R')
source('R/diversity.r')
source('R/correlationTests.r')
source('R/taxaplot.R')

# ----- READ DATA AND CALCULATE BASIC STATS -----

# Read metadata
md = read.csv(file="metadata/metadata.csv",header=T,row.names=1)
row.names(md) = gsub("-",".",row.names(md))
mdo= md[order(row.names(md)),]
summary(mdo)

# Type     Nacid     Rep         Run            Primer     
# ctr:  4   DNA:439   1:448   2017  : 51   515..926 : 51  
# S  :260   RNA: 29   A:  10   201809:  9   519..806 : 27  
# W  :204             B: 10    201812: 190  519..806D : 390
#                             201904: 191
#                             Eva: 27
# 
# 375 est, 46 offshore, 39 coastal and 4 port. AMBI .5 -- 5.991

summary(as.factor(mdo$Year))
# 2013 2015 2016 2017 2018      2019   NA
# 64    25    60   136   125      54    4

otu.file = paste(otudir, "SWARM_table_curated.tsv", sep="/")
otus.all = read.delim(otu.file,row.names=1,header=T,sep="\t")
dim(otus.all)
# 188224, 469

otus.all = otus.all[,order(names(otus.all))]
classification=otus.all$classification
otus.t = as.data.frame(t(otus.all[,-3]))

table(row.names(mdo) %in% rownames(otus.t))
# and do they match exactly?
table(row.names(otus.t) == row.names(mdo))

write.csv(rowSums(otus.t),"Diversity_stats/reads_in_OTUs.csv")

sum(otus.t)
#16,045,227 reads
summary(rowSums(otus.t))
# Min.  125 

hist(rowSums(otus.t),breaks=100,xlim=range(0,50000))

otus.clean = filterCrossContaminants2(otus.t,100)
sum(otus.clean)
# 16,041,964

write.csv(rowSums(otus.t)-rowSums(otus.clean),"Diversity_stats/CrossContaminantReads.csv",quote=F)

# Remove unclassified 
unclassified = otus.all$classification=="No hits"
otus.class = otus.clean[,!unclassified]
classified = otus.all$classification[!unclassified]
names(classified) = names(otus.class)

sum(otus.class) #15,994,075 (99.7%)
dim(otus.class) #185,694 OTUs (98.7%)

# ------- REMOVE CONTAMINANTS BASED ON CONTROL ------

type = paste(mdo$Type, mdo$Run, mdo$Environment, mdo$Nacid)
nType = length(unique(type))
pooledByType = matrix(nrow = nType, ncol = dim(otus.class)[2])
for (i in 1:nType){
  pool = unique(type)[i]
  pooledByType[i,] = colSums(otus.class[type==pool,])
}
pooledByType = data.frame(pooledByType)
rownames(pooledByType) = unique(type)
names(pooledByType) = names(otus.class)

rowSums(pooledByType)

write.table(t(pooledByType), "pooledByType.tsv", sep="\t",quote=F)

otus.blank = otus.class[mdo$Type=="ctr",]
otus.temp = rbind(otus.blank,pooledByType)
otus.ra.temp = decostand(otus.temp,method="total")

otus.inBlank = otus.ra.temp[,colSums(otus.blank) > 0]

dim(otus.inBlank)
# 509 OTUs are present in blanks

pdf("img/contamination/OTUs_in_controls_barchart.pdf",height=12,width=14)
taxaplot(25,data.frame(row.names=row.names(otus.inBlank),
                       c(rep("Controls",5),rep("Samples",17))),otus.inBlank)
dev.off()

# Clear contamination candidates: SWARMs:
# Extraction: 2164 5912 744
# PCR: 7? 227 119 8 53

classified["SWARM_2164"] # Vibrio
classified["SWARM_5912"] # Streptococcus
classified["SWARM_744"] # Acinetobacter
classified["SWARM_227"] # Vibrio

library(matrixStats)
contRatio = t(otus.inBlank[5,] / colMaxs(as.matrix(otus.inBlank[c(6:22),])))
hist(log10(contRatio),breaks=50,col="grey",
     main="Abundance in control / max. abundance by type")

trueContaminants = otus.inBlank[,contRatio>10]
unContaminants = otus.inBlank[,contRatio<=10]
dim(trueContaminants) #200 of 510

pdf("img/contamination/real_contaminants_barchart.pdf",height=12,width=14)
taxaplot(15,data.frame(row.names=row.names(trueContaminants),
                       c(rep("Controls",5),rep("Samples",17))),trueContaminants)
dev.off()

pdf("img/contamination/OTUs_in_blank_retained.pdf",height=12,width=14)
taxaplot(15,data.frame(row.names=row.names(unContaminants),
                       c(rep("Controls",5),rep("Samples",17))),unContaminants)
dev.off()
contaminant_list = names(trueContaminants)
pooledByType.clean = pooledByType[,!(names(pooledByType) %in% contaminant_list)]
dim(pooledByType.clean)

classified.clean = classified[!(names(classified) %in% contaminant_list)]
otus.class.clean = otus.class[,!(names(otus.class) %in% contaminant_list)]

writeDivStats("Diversity_stats/Diversity_by_type.csv",pooledByType.clean[-1,])

# -------- EUKARYOTES - Remove and save separately ----------------

euk = grep("Eukaryota",classified.clean)
otus.p = otus.class.clean[,-euk]
classified.p = classified.clean[-euk]
sum(otus.p) #15,402,002 (96.5%)
sum(otus.p)/sum(otus.class.clean)
dim(otus.p) #182,705 OTUs

otus.e = otus.class.clean[,euk]
sum(otus.e) #550,864 reads
dim(otus.e) #2789 OTUs

chloro = grep("Chloroplast;",classified.clean)
chloroTax = classified.clean[chloro]
otus.chl = otus.class.clean[,chloro]
sum(otus.chl) #526,395s
dim(otus.chl) #1116 

write.table(t(otus),paste(otudir, "../Chloroplast_OTU_table_final.tsv", sep="/"), sep="\t",quote=F,
            col.names=NA)

# Most abundant unclassified
otus.uc = otus.t[,unclassified]
otus.uc = otus.uc[,order(colSums(otus.uc),decreasing = T)]
names(otus.uc)[c(1:10)]
colSums(otus.uc)[c(1:10)]

# Most abundant euk
otus.e = otus.e[,order(colSums(otus.e),decreasing = T)]
names(otus.e)[c(1:10)]
classified.clean[names(otus.e)[c(1:10)]]
colSums(otus.e)[c(1:10)]

# QC statistics
write.csv(rowSums(otus.clean)-rowSums(otus.class),"Diversity_stats/UnclassifiedReads.csv",quote=F)
write.csv(rowSums(otus.class)-rowSums(otus.class.clean),"Diversity_stats/ContaminantReads.csv",quote=F)
write.csv(rowSums(otus.class.clean)-rowSums(otus.p),"Diversity_stats/EukReads.csv",quote=F)

sum(rowSums(otus.p)>5000) #456 / 486
sum(rowSums(otus.p)>10000) #425
otus = otus.p[rowSums(otus.p)>10000,]
mdi = mdo[rowSums(otus.p)>10000,]

# Remove LN10 that was misidentified!
otus=otus[-c(row.names(mdi)=="LN10_S_Wi2013_PS"),]
mdi=mdi[-c(row.names(mdi)=="LN10_S_Wi2013_PS"),]

otus = otus[,colSums(otus)>0] 
dim(otus)# 182502
taxonomy = classified.p[colSums(otus)>0]

# Div statistics

writeDivStats("Diversity_stats/Diversity_classified_all.csv", otus.p)
writeDivStats("Diversity_stats/Diversity_classified.csv", otus)
writeDivStats("Diversity_stats/Diversity_euk.csv", otus.e)
divTemp = read.csv("Diversity_stats/Diversity_classified.csv", row.names=1, header=T)

lmrs=lm(Richness~Reads,data=divTemp)
summary(lmrs) # p=1E-14
pdf("img/Richness_vs_reads.pdf")
plot(Richness~Reads,data=divTemp,col=as.numeric(mdi$Type)+2)
abline(lmrs,col="grey",lty=2)
legend("bottomright",col=c(4,3),legend=c("Water","Sediment"),pch=1)
dev.off()

lmrrs=lm(Rarefied.richness~Reads,data=divTemp)
summary(lmrrs) # p=0.4
pdf("img/Rarefied_richness_vs_reads.pdf")
plot(Rarefied.richness~Reads,data=divTemp,col=as.numeric(mdi$Type)+2)
abline(lmrrs,col="grey",lty=2)
legend("bottomright",col=c(4,3),legend=c("Water","Sediment"),pch=1)
dev.off()

# --------- Write cleaned-up OTU and metadata tables -------

write.table(t(otus),paste(otudir, "../SWARM_table_final.tsv", sep="/"), sep="\t",quote=F,
            col.names=NA)

write.table(mdi,"metadata/metadata_final.tsv", sep="\t",quote=F, col.names=NA)

