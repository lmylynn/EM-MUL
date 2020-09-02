# EM-MUL
EM-MUL is an effective tools which resolves ambiguous bisulfite-treated reads, making use of information we have.  
#To run this program, we needs to have samtools, perl, bedtools,and g++ first.  
The  inputs of this tool consists four parts. 
-r is the reference genome to be aligned.
-u is the unique reads.
-m is multireads,which align to multiple locations of the reference genome ambiguously.
-o is the unique reads that overlapped with multireads.
Among them, the unique reads and the multireads are obtained by aligning the original BS reads to bismark. 
Overlappedfile can be obtained through the unique reads and multireads, the processing flow refers to BAM_ABS, the commad is: 
1. Convert unique_reads.sam to unique_reads.bam.
  samtools view -bS unique_reads.sam > unique_reads.bam 
2. Run Covert_to_bed_unite.pl to covert ambiguous read file to bed formate with --ambiguous option.
  perl Convert_to_bed_unite.pl --ambiguous ambiguous_file.sam 
3.Run samtools to get overlapped unique reads in sam format. 
  samtools view -L ambiguous_file.bed unique_reads.bam -q 20 > unique_reads.sam 
4.To get rid of duplicates from the unique reads.
  sort -n -r -k3,3 -k4,4 -k5,5 unique_reads.sam|uniq -u > unique_reads_nodup.sam 
5. Convert unique read file to bed format with --unique option.
  perl Convert_to_bed_unite.pl --unique unique_reads_nodup.sam 
6. Get overlapped unique reads by using Bedtools and run the following command in the bedtools folder to get the overlappedfile we use.
  ./intersectBed -a ambiguous_file.bed -b unique_reads_nodup.bed -wb -wa > overlapfile.txt 
7. Score the multireads using EM-MUL.
  #python3 new_score_all_and_coverage_human -r hg38 -u unique_reads.sam -m multireads.sam -r overlapfile.sam
