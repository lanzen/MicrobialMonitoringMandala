setwd(".")

require(lulu)
require(methods)

matchlist = read.delim("LULU_match_list.txt", header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
otus.all = read.delim("CREST_Results/SWARM_table.tsv",row.names=1,header=T,sep="\t")
otus = otus.all[,-dim(otus.all)[2]]
curated_result <- lulu(otus,matchlist, minimum_match = 97)
lulus = curated_result$curated_table
write.table(data.frame("OTU"=rownames(lulus),lulus),"SWARM_table_curated.tsv", row.names=FALSE, 
            quote=F, sep="\t")
