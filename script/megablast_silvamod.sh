blastn -task megablast -query $1 -db ~/projects/PhyloRefDB/silvamod128/silvamod128.fasta -num_alignments 100 -outfmt 5 -out ${1//.fa*}_silvamod128.xml -num_threads 3
gzip ${1//.fa*}_silvamod128.xml
