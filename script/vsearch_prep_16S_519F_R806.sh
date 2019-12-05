# v 43 2017-09-08 mod 2018-10-03 (RD and liberal criteria)

fprimer=CAGCMGCCGCGGTAA
fprimlength=15
#C:  GGACTAC--H--VGGGTWTCTAAT
#RD: GGACTAC->N<-VGGGTWTCTAAT
#
rprimer_rc=ATTAGAWACCCBNGTAGTCC

#Length of amplicon minus primers (787-519-15= 253)
#(the rev. primer is 20 nt but starts at pos. 787)
crop_length=252 

#Make directory where the merged read files will be placed (fastq and fasta)
mkdir merged_reads

#Iterate over all subdirectories, assume they have the sample name, enter and
#uncompress fastq files
for d in $*; do
    echo $d
    cd $d
    pwd
    i=${d//\/}
    echo $i

    gunzip -kf *.fastq.gz 
    
    # Merge the reads (non-staggering since adaptors removed) max 20 differences (defualt 5)

    echo $d >> ../readprep.log

    vsearch --fastq_mergepairs *_R1_001.fastq --reverse *_R2_001.fastq  --fastq_allowmergestagger --fastqout ${i}_m.fastq --fastq_maxdiffs 20 2>> ../readprep.log

    # Allow up to 2 mismatches and one missing base in F primer. still only about 33%
    cutadapt -g $fprimer -m 100 --max-n=0 --discard-untrimmed -O 14 -e 0.15 -o ${i}_trimmed1.fastq  ${i}_m.fastq >> ../readprep.log 

    cutadapt -a $rprimer_rc -m 100 --discard-untrimmed -O 19 -o ${i}_trimmed.fastq ${i}_trimmed1.fastq >>../readprep.log

    rm  ${i}_trimmed1.fastq

    # Max expected errors 1 (default 0.5)
    vsearch --fastq_filter ${i}_trimmed.fastq --fastaout ${i}_merged_QF.fasta --fastq_maxee 1 --fastq_trunclen $crop_length 2>>../readprep.log    
    
    #Fix read names to include "barcode" label
    fixreads.py ${i}_merged_QF.fasta $i > ../merged_reads/${i}_reads_fixed.fasta

   rm ${i}_merged_QF.fasta

   cd ..
done

echo "Preparing FASTQC reports"

cat */*_m.fastq > all_merged.fastq
fastqc all_merged.fastq
rm all_merged.fastq

cat */*_trimmed.fastq > all_trimmed.fastq
fastqc all_trimmed.fastq
rm all_trimmed.fastq

cat */*R1*.fastq > all_R1.fastq
fastqc all_R1.fastq
rm all_R1.fastq

cat */*R2*.fastq > all_R2.fastq
fastqc all_R2.fastq
rm all_R2.fastq

rm */*.fastq
rm *fastqc.zip

for h in *.html; do gnome-open $h; done

