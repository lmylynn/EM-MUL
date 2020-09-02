EM-MUL
====
EM-MUL is an effective tools which resolves ambiguous bisulfite-treated reads, making use of information we have.
To run this program, we needs to have samtools, perl, bedtools,and g++ first.  <br>
The  inputs of this tool consists four parts. <br>
* -r is the reference genome to be aligned.<br>
* -u is the unique reads.<br>
* -m is multireads,which align to multiple locations of the reference genome ambiguously.<br>
* -o is the unique reads that overlapped with multireads.<br>

Among them, the unique reads and the multireads are obtained by aligning the original BS reads to bismark. 
Overlappedfile can be obtained through the unique reads and multireads, the processing flow refers to BAM_ABS, the commad is: <br>
* Convert unique_reads.sam to unique_reads.bam.<br>
    * samtools view -bS unique_reads.sam > unique_reads.bam <br> 
* Run Covert_to_bed_unite.pl to covert ambiguous read file to bed formate with --ambiguous option.<br>
    * perl Convert_to_bed_unite.pl --ambiguous ambiguous_file.sam <br>
* Run samtools to get overlapped unique reads in sam format. <br>
    * samtools view -L ambiguous_file.bed unique_reads.bam -q 20 > unique_reads.sam <br>
* To get rid of duplicates from the unique reads.<br>
    * sort -n -r -k3,3 -k4,4 -k5,5 unique_reads.sam|uniq -u > unique_reads_nodup.sam <br> 
* Convert unique read file to bed format with --unique option.<br>
    * perl Convert_to_bed_unite.pl --unique unique_reads_nodup.sam <br>
* Get overlapped unique reads by using Bedtools and run the following command in the bedtools folder to get the overlappedfile we use.<br>
    * ./intersectBed -a ambiguous_file.bed -b unique_reads_nodup.bed -wb -wa > overlapfile.txt <br>
* Score the multireads using EM-MUL.<br>
    * python3 new_score_all_and_coverage_human -r hg38 -u unique_reads.sam -m multireads.sam -r overlapfile.sam<br>
