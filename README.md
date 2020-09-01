# EM-MUL
EM-MUL is an effective tools which resolves ambiguous bisulfite-treated reads, making use of information we have.  
#To run this program, we needs to have samtools, perl, bedtools,and g++ first.  
The  inputs of this tool consists four parts. 
-r is the reference genome to be aligned.
-u is the unique reads.
-m is multireads,which align to multiple locations of the reference genome ambiguously.
-o is the unique reads that overlapped with multireads.
Among them, the unique reads and the multireads are obtained by aligning the original BS reads to bismark. 
Overlapfile can be obtained through the unique reads and multireads, the processing flow refers to BAM_ABS, the commad is: 
samtools view -bS unique_reads.sam > unique_reads.bam 
perl Convert_to_bed_unite.pl --ambiguous ambiguous_file.sam 
samtools view -L ambiguous_file.bed unique_reads.bam -q 20 > unique_reads.sam 
sort -n -r -k3,3 -k4,4 -k5,5 unique_reads.sam|uniq -u > unique_reads_nodup.sam 
perl Convert_to_bed_unite.pl --unique unique_reads_nodup.sam 
./intersectBed -a ambiguous_file.bed -b unique_reads_nodup.bed -wb -wa > overlapfile.txt 
#python3 new_score_all_and_coverage_human -r hg38 -u unique_reads.sam -m multireads.sam -r overlapfile.sam
