# ----- READ DATA -----

setwd("~/projects/uGAMBI_RC/MBC_2019/")
otudir = "SWARM_20190503/CREST_LULU_Final"
require(vegan)

source('R/filtering.R')
source('R/diversity.r')
source('R/correlationTests.r')
source('R/taxaplot.R')

mdi = read.table(file="metadata/metadata_final.tsv",header=T,row.names=1,sep="\t")
summary(mdi)
mdf = mdi[,c("Type","SalClass")] #factors
mdc = mdi[,c("Year","microgAMBI","PI","MaxDepthPM","PointDepthBM","Salinity",
             "SalClass","Q1SiteSal","Q2SiteSal","Q3SiteSal")]

otu.file = paste(otudir, "SWARM_table_final.tsv.gz", sep="/")
otus.t = read.delim(otu.file,row.names=1,header=T,sep="\t")
dim(otus.t)

otus.t = otus.t[,order(names(otus.t))]
taxonmoy=otus.t$classification
otus = as.data.frame(t(otus.t[,-1]))

# Are all samples represented?
table(row.names(mdi) %in% rownames(otus))
# and do they match exactly?
table(row.names(otus) == row.names(mdi))

divTemp = read.csv("Diversity_stats/Diversity_classified.csv", row.names=1, header=T)
div16S = divTemp[,c(5:7,9:10)]
summary(divTemp)
# Ranging from 10,141 reads to 138,147. Richness 290 -- 13,384, rarefied: 235 -- 4623.. 
row.names(div16S) == row.names(mdi)

# ---- Calculate relative abundance and Hellinger transformed and filter rare ----
# ( sequencing depth uneven with smallest sample having 10000 reads )

otus.ra.all = decostand(otus, method="total")
otus.ra.f.16S1 = dropRareByMaxAbundance(otus.ra.all, 10/min(rowSums(otus)))
otus.ra.f.16S = decostand(otus.ra.f.16S1, method="total")
dim(otus.ra.f.16S) # 5855 OTUs 
otus.H = decostand(otus.ra.f.16S1, method="hell")
otus.ra.bc = vegdist(otus.ra.f.16S)

# Graph colours
rb = rainbow(17)
colour = rb[as.numeric(mdi$Estuary)] 
colU = rb[sort(as.numeric(unique(mdi$Estuary)))]
# make None = black
colour[mdi$Estuary == "None"] = "black"
colU[sort(unique(mdi$Estuary)) == "None"] = "black"

# ------- Alphadiv correlation -----

printANOVA(mdf, div16S, a=0.05)
# Diversity of sediments clearly higher than water by all measures

printVS(div16S,mdc, a=.05/10)
# Rarefied richness, H and J neg. corr w microgAMBI across both types

lmmgs=lm(div16S$Rarefied.richness~mdi$microgAMBI)
summary(lmmgs) # p=1E-11, R2=.10,  t=-6.8 (-398+-58)
plot(div16S$Rarefied.richness~mdi$microgAMBI,col=as.numeric(mdi$Type)+2,xlab="microgAMBI",
     ylab="Rarefied richness")
# The trend is in fact opposed with water of high microgAMBI having higher div!

lmmgj=lm(div16S$J~mdi$microgAMBI)
summary(lmmgj) # p=8E-14, R2=.12,  t=-7.7
pdf("img/Evenness_vs_microgAMBI.pdf")
plot(div16S$J~mdi$microgAMBI,col=as.numeric(mdi$Type)+2,xlab="microgAMBI",
     ylab="Pielou evenness")
abline(lmmgj,col="grey",lty=2)
legend("bottomright",col=c(4,3),legend=c("Water","Sediment"),pch=1)
dev.off()


# ------- ANOSIM and adonis --------

anosim(otus.ra.bc, grouping=mdi$Type)
# (water vs. sediment) R = 0.88, p<0.001
anosim(otus.H, grouping=mdi$Type)
# R=0.90, slightly more sensitive

adonis(otus.H~Type+Environment+Season+Station+Year+Primer,data=mdi)

# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Type          1    38.053  38.053 265.676 0.27044  0.001 ***
#   Environment   3     9.236   3.079  21.495 0.06564  0.001 ***
#   Season        3     6.948   2.316  16.170 0.04938  0.001 ***
#   Station      67    35.266   0.526   3.675 0.25064  0.001 ***
#   Year          1     0.413   0.413   2.882 0.00293  0.009 ** 
#   Primer        2     1.088   0.544   3.800 0.00774  0.001 ***
#   Residuals   347    49.701   0.143         0.35323           
# Total       424   140.706                 1.00000           
         


# -------- NMDS --------

# Hellinger transformed has also been tried here but makes very little difference

nmdsRA = metaMDS(otus.ra.f.16S) 
eff = envfit(nmdsRA, mdc,na.rm=T)
eff

 
# NMDS1    NMDS2     r2 Pr(>r)    
# Year         -0.99504  0.09943 0.1327  0.001 ***
#   microgAMBI   -0.18779  0.98221 0.5932  0.001 ***
#   PI            0.99877  0.04965 0.6099  0.001 ***
#   MaxDepthPM   -0.24678 -0.96907 0.1662  0.001 ***
#   PointDepthBM -0.35000 -0.93675 0.1714  0.001 ***
#   Salinity     -0.03020 -0.99954 0.3786  0.001 ***
#   SalClass     -0.03758 -0.99929 0.4036  0.001 ***
#   Q1SiteSal    -0.00439 -0.99999 0.2065  0.001 ***
#   Q2SiteSal    -0.12053 -0.99271 0.2090  0.001 ***
#   Q3SiteSal    -0.26992 -0.96288 0.1940  0.001 ***
#   ---

effF = envfit(nmdsRA, mdf,na.rm=T)
#   Type    0.6525  0.001 ***
#   Estuary 0.2026  0.001 ***
#   Tide    0.2818  0.001 ***
#   Season  0.3188  0.001 ***

pdf("img/NMDS_type_coloured.pdf")
ordiplot(nmdsRA,type="none")
points(nmdsRA,display="sites", pch=as.numeric(mdi$Type),lwd=2,cex=.7,
       col=as.numeric(mdi$Type)+2)
legend("bottomleft",col=c(4,3),legend=c("Water","Sediment"),pch=c(2,1),cex=.7)
dev.off()

pdf("img/NMDS_est_coloured.pdf")
ordiplot(nmdsRA,type="none")
points(nmdsRA,display="sites",col=colour, pch=as.numeric(mdi$Type),lwd=2,cex=.8)
legend("bottomleft",col=colU, legend=sort((unique(mdi$Estuary))),pch=1,cex=.7,ncol=5)
dev.off()

# ------ TAXONOMIC BARCHARTS --------

# Read taxonomic data
taxa.all = read.delim(paste(otudir, "Relative_Abundance.tsv.gz", sep="/"), header=T, row.names=3)
table(row.names(mdi)==names(taxa.all)[-c(1:2)]) # yes
t
#grouping_info<-data.frame(row.names=names(taxa.all)[-c(1:2)], mdi$Estuary)
grouping_info<-data.frame(row.names=row.names(mdi), mdi$Type)

ranks = data.frame(rank=c("domain","superkingdom","kingdom","phylum","class","order","family","genus","species"),
                   levels=c(3,10,21,21,21,21,21,21,21))

for (i in c(1:dim(ranks)[1])){
  level=toString(ranks$rank[i])
  levelTaxa = as.data.frame(t(taxa.all[taxa.all$Rank==level,-c(1:2)]))
  print(paste(mean(rowSums(levelTaxa))*100,"% classified at rank",level))
  #line <- readline()      
  pdf(paste("img/16S_barcharts/",level,"_barchart.pdf",sep=""),height=10,width=26)
  taxaplot(ranks$levels[i],grouping_info,levelTaxa)
  dev.off()
}


# ------ Distribution of OTUs across environs -------

otus.s = otus[mdi$Type=="S",]
otus.s = otus.s[,colSums(otus.s)>0]

otus.w = otus[mdi$Type=="W",]
otus.w = otus.w[,colSums(otus.w)>0]

dim(otus.s)[2]+dim(otus.w)[2]-dim(otus)[2]
# 39,794 shared between both environment types (22%)

dim(otus)[2]-dim(otus.w)[2]
# 124,584 unique to sed
dim(otus)[2]-dim(otus.s)[2]
# 18,127 unique to water


