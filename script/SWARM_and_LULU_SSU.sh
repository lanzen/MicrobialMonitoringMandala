#!/bin/bash

# Anders Lanzen 2019-05-03 v1.1
# This script uses vsearch for dereplication of each individual sample
# then merges them using code from SLIM, retaining read mapping, allowing
# for an OTU table to be constructed after SWARM clustering. Singletons are
# removed by default (change minAbundance to 1 to retain) and chimeras removed
# using Silvamod and de novo.  Rmoved OTUs are
# removed from the OTU table.

minAbundance=2 # Change to 1 to retain singletons
db_dir=$HOME/projects/PhyloRefDB/silvamod128
scriptdir=$HOME/script


#Dereplicate each file separately instead of concatenating at first
for fa in *_fixed.fa*; do 
    sed 's/ 1:N:0:.*;/;/g' $fa > ${fa//.*}_f.fasta
    vsearch -derep_fulllength ${fa//.*}_f.fasta -output ${fa//_reads_fixed.fa*}.dfa -sizeout
    rm ${fa//.*}_f.fasta
done
 
fasta_merging.py *.dfa
gzip *.dfa

#Abundance sort and retain singletons 
vsearch -sortbysize merged.fasta -output merged_sorted.fasta -minsize 0
# rm merged.fasta

# Run SWARM
swarm -f -z -t 4 -w SWARM_OTUs.fasta merged_sorted.fasta -o SWARM.swarms -s SWARM.stats -u SWARM.uc

# Create OTU table

matrix_creation.py -uc SWARM.uc -so origins.tsv -o SWARM_table.tsv -fasta_in merged_sorted.fasta -fasta_out merged_swarm.fasta
rm merged_sorted.fasta

# Annotate SWARM table and OTUs with number and "SWARM_" prefix
sed -i -e 's/^/SWARM_/' SWARM_table.tsv 
fasta_number.py SWARM_OTUs.fasta SWARM_ -needsize > SWARM_OTUs_f.fasta

# Sort by size and remove singleton swarms
vsearch --sortbysize SWARM_OTUs_f.fasta -output SWARM_OTUs_ns.fasta -minsize 2

# Chimera filtering
vsearch --uchime_ref SWARM_OTUs_ns.fasta -db $db_dir/silvamod128.fasta -strand plus -nonchimeras SWARM_OTUs_uchime_ref.fa

vsearch --uchime_denovo SWARM_OTUs_uchime_ref.fa --nonchimeras SWARM_OTUs_uchime_ref_denovo.fa

awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' SWARM_OTUs_uchime_ref_denovo.fa  > SWARM_OTUs_L.fasta

sed 's/;size=.*//' SWARM_OTUs_L.fasta > SWARM_OTUs_clean.fasta
rm SWARM_OTUs_L.fasta

# Classify using CREST which also removes discarded OTUs

filter_OTUtable_by_seq.py SWARM_table.tsv SWARM_OTUs_clean.fasta > SWARM_table_clean.tsv


# Align for LULU post-clustering

makeblastdb -in SWARM_OTUs_clean.fasta -dbtype nucl

blastn -db SWARM_OTUs_clean.fasta -outfmt '6 qseqid sseqid pident' -out LULU_match_list.txt -qcov_hsp_perc 95 -perc_identity 95 -query SWARM_OTUs_clean.fasta -num_threads 3

Rscript ~/script/LULU.R

# Classify LULU filtered OTUs using CREST


filter_seq_by_OTUtable.py SWARM_OTUs_clean.fasta SWARM_table_curated.tsv > SWARM_OTUs_curated.fasta

megablast_silvamod.sh SWARM_OTUs_curated.fasta

classify -o CREST_LULU -i SWARM_OTUs_curated.fasta -t SWARM_table_curated.tsv SWARM_OTUs_curated_silvamod128.xml.gz



