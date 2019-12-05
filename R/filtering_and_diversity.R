setwd("~/projects/uGAMBI_RC/Mandala_Sediments_Article1") #<- change to correct local directory

otudir = "SWARM_20190911/CREST_LULU"

require(vegan)
source('R/utils/filtering.R')
source('R/utils/diversity.r')
source('R/utils/correlationTests.r')
source('R/utils/taxaplot.R')

# ----- READ DATA AND CALCULATE BASIC STATS -----

# Read metadata and fix names
md = read.csv(file="metadata/metadata.csv",header=T,row.names=1)
row.names(md) = gsub("-",".",row.names(md))
mdo= md[order(row.names(md)),]

# Read OTU table
otu.file = paste(otudir, "SWARM_table_curated.tsv.gz", sep="/")
otus.all = read.delim(otu.file,row.names=1,header=T,sep="\t")
dim(otus.all)
# 215,755  600

# Order OTUs by name, transverse (for vegan) and control correspodnence to metadata
otus.all = otus.all[,order(names(otus.all))]
classification=otus.all$classification
otus.t = as.data.frame(t(otus.all[,-3]))

data.frame(row.names(otus.t),row.names(mdo))[row.names(otus.t) != row.names(mdo),]
table(row.names(otus.t) == row.names(mdo))

# ------ Basic reads stats -------
write.csv(rowSums(otus.t),"Diversity_stats/reads_in_OTUs.csv")

sum(otus.t)
#23,031,579 reads
summary(rowSums(otus.t))
# Min.  54, Q1 23090 

hist(rowSums(otus.t),breaks=100,xlim=range(0,50000))

# ------ Cross contamination and unclassified filtering -------

otus.clean = filterCrossContaminants2(otus.t,100)
sum(otus.clean)
# 23,027,227

write.csv(rowSums(otus.t)-rowSums(otus.clean),"Diversity_stats/CrossContaminantReads.csv",quote=F)

# Remove unclassified 
unclassified = otus.all$classification=="No hits"
otus.class = otus.clean[,!unclassified]
classified = otus.all$classification[!unclassified]
names(classified) = names(otus.class)

sum(otus.class) #22,952,976 (99.%)
dim(otus.class) #212,688 OTUs (98.6%)


# Most abundant unclassified
otus.uc = otus.t[,unclassified]
otus.uc = otus.uc[,order(colSums(otus.uc),decreasing = T)]
names(otus.uc)[c(1:10)]
colSums(otus.uc)[c(1:10)]

#SWARM_202 and 550: probable Ostreococcus spp. mitochondrion 
#SWARM_1873: Bathycoccus prasinos mito 
# (both have chloroplast seqs classified)

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
# 2202 OTUs are present in blanks

pdf("img/contamination/OTUs_in_controls_barchart.pdf",height=12,width=14)
taxaplot(25,data.frame(row.names=row.names(otus.inBlank),
                       c(rep("Controls",11),rep("Samples",23))),otus.inBlank)
dev.off()

# Clear contamination candidates: SWARMs:
# Extraction: 2164 5912 744
# PCR: 7? 227 119 8 53

classified["SWARM_2807"] # Vibrio
classified["SWARM_4818"] # Streptococcus
classified["SWARM_382"] # Chloroplast
classified["SWARM_30"] # Pseudoalteromonas

library(matrixStats)
contRatio = t(colMaxs(as.matrix(otus.inBlank[c(1:8),]))
              / colMaxs(as.matrix(otus.inBlank[c(12:34),])))
hist(log10(contRatio),breaks=50,col="grey",
     main="Max. abundance in control / max. abundance by type")

trueContaminants = otus.inBlank[,contRatio>10]
unContaminants = otus.inBlank[,contRatio<=10]
dim(trueContaminants) #1006 of 2022

pdf("img/contamination/real_contaminants_barchart.pdf",height=12,width=14)
taxaplot(15,data.frame(row.names=row.names(trueContaminants),
                       c(rep("Controls",11),rep("Samples",23))),trueContaminants)
dev.off()

pdf("img/contamination/OTUs_in_blank_retained.pdf",height=12,width=14)
taxaplot(15,data.frame(row.names=row.names(unContaminants),
                       c(rep("Controls",11),rep("Samples",23))),unContaminants)
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
sum(otus.p) #21,543,065 (95%)
sum(otus.p)/sum(otus.class.clean)
dim(otus.p) #207,897  OTUs

otus.e = otus.class.clean[,euk]
sum(otus.e) #1,120,894 reads
dim(otus.e) #3785 OTUs

chloro = grep("Chloroplast;",classified.clean)
chloroTax = classified.clean[chloro]
otus.chl = otus.class.clean[,chloro]
sum(otus.chl) #1,083,344
dim(otus.chl) #1497 

# Write OTU data
write.table(t(otus.e),paste(otudir, "../Euk_OTU_table_final.tsv", sep="/"), sep="\t",quote=F,
            col.names=NA)

# Water only for chloroplasts
otus.w.chl = otus.chl[mdo$Type=="W",]
otus.w.chl = otus.w.chl[,colSums(otus.w.chl)>0]
otus.w.p = otus.p[mdo$Type=="W",]
otus.w.p = otus.w.p[,colSums(otus.w.p)>0]
mdw = mdo[mdo$Type=="W",]
write.table(mdw,"metadata/metadata_water_all.tsv", sep="\t",quote=F, col.names=NA)

write.table(t(otus.w.chl),paste(otudir, "../Chloro_OTU_table_water.tsv", sep="/"), sep="\t",quote=F,
            col.names=NA)

write.table(t(otus.w.p),paste(otudir, "../OTU_table_water_all.tsv", sep="/"), sep="\t",quote=F,
            col.names=NA)

# Most abundant euk
otus.e = otus.e[,order(colSums(otus.e),decreasing = T)]
names(otus.e)[c(1:10)]
classified.clean[names(otus.e)[c(1:10)]]
colSums(otus.e)[c(1:10)]
# Only classified as "chloroplast" period


# QC statistics
write.csv(rowSums(otus.clean)-rowSums(otus.class),"Diversity_stats/UnclassifiedReads.csv",quote=F)
write.csv(rowSums(otus.class)-rowSums(otus.class.clean),"Diversity_stats/ContaminantReads.csv",quote=F)
write.csv(rowSums(otus.class.clean)-rowSums(otus.p),"Diversity_stats/EukReads.csv",quote=F)

# ------ Remove samples with too low prok. sequencing depth -------

sum(rowSums(otus.p)>5000) #582 / 599
sum(rowSums(otus.p)>10000) #569 / 599
otus = otus.p[rowSums(otus.p)>10000,]
rowSums(otus.chl)
mdi = mdo[rowSums(otus.p)>10000,]

# Remove LN10 that was misidentified!
otus=otus[c(row.names(mdi)!="LN10_S_Wi2013_PS"),]
mdi=mdi[c(row.names(mdi)!="LN10_S_Wi2013_PS"),]

otus = otus[,colSums(otus)>0] 
dim(otus)# 207,724
taxonomy = classified.p[colSums(otus)>0]

# ------ Div statistics -------

writeDivStats("Diversity_stats/Diversity_classified.csv", otus)
writeDivStats("Diversity_stats/Diversity_euk.csv", otus.e)
divTemp = read.csv("Diversity_stats/Diversity_classified.csv", row.names=1, header=T)
divp = divTemp[,c("Reads","Richness","Rarefied.richness","H","J","Chao1")]

lmrs=lm(Richness~Reads,data=divp)
summary(lmrs) # p=7E-3
pdf("img/Richness_vs_reads.pdf",width=6,height=6)
plot(Richness~Reads,data=divTemp,col=as.numeric(mdi$Type)+1,
     pch=as.numeric(mdi$Nacid))
abline(lmrs,col="grey",lty=2)
legend("bottomright",col=c(4,3,1,1),legend=c("Water","Sediment","DNA","RNA"),
       pch=c(1,1,1,2))
dev.off()

lmrrs=lm(Rarefied.richness~Reads,data=divTemp)
summary(lmrrs) # p=6E-7, t=-5 but appears to be effect of more sequenced water samples
pdf("img/Rarefied_richness_vs_reads.pdf")
plot(Rarefied.richness~Reads,data=divTemp,col=as.numeric(mdi$Type)+2)
abline(lmrrs,col="grey",lty=2)
legend("bottomright",col=c(4,3),legend=c("Water","Sediment"),pch=1)
dev.off()

# --------- Write cleaned-up OTU and metadata tables -------

write.table(t(otus),paste(otudir, "../SWARM_table_final.tsv", sep="/"), sep="\t",quote=F,
            col.names=NA)

write.table(mdi,"metadata/metadata_final.tsv", sep="\t",quote=F, col.names=NA)


