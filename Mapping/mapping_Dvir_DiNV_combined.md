# Mapping Reads to a Concatenated _D. virilis_ and DiNV Genome

**Step 1: Create Concatenated Genome**
- Right now I have the _D. virilis_ genome as multiple chromosomes, and I have the DiNV genome separate as well
- To map my sequences to these, I'll have to combine/concatenate them   
`cat *fna > Dvir.DiNV.combo.fa`
- Now I want to make a new directory where I can do the mapping  
`mkdir Mapping`  
- Next I want to symbolically link (or symlink) my concatenated genome to the new folder so I don't have to copy it and use up more space  
`ln -s /home/runcklesslab/Maggie/Dvir.DiNV.combo.fa /home/runcklesslab/Maggie/Mapping/Dvir.DiNV.combo.fna`
- And I need the trimmed sequence files in this directory too  
`ln -s /home/runcklesslab/Maggie/KM_3_1_trim_d.fq.gz Mapping/KM_1_trim_d.fq.gz`  
`ln -s /home/runcklesslab/Maggie/KM_3_2_trim_d.fq.gz Mapping/KM_2_trim_d.fq.gz`

**Step 2: Using BWA to Create a Genome Index and Map to the Genome**
- [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) is an alignment tool that is widely used
- First "BWA first needs to construct the FM-index for the reference genome"  
`bwa index Dvir.DiNV.combo.fna`  
- Output from bwa index:   

```
[bwa_index] Pack FASTA... 1.03 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=339854600, availableWord=35912944
[BWTIncConstructFromPacked] 10 iterations done. 59240328 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 109442536 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 154058360 characters processed.
[BWTIncConstructFromPacked] 40 iterations done. 193708968 characters processed.
[BWTIncConstructFromPacked] 50 iterations done. 228946536 characters processed.
[BWTIncConstructFromPacked] 60 iterations done. 260261704 characters processed.
[BWTIncConstructFromPacked] 70 iterations done. 288090696 characters processed.
[BWTIncConstructFromPacked] 80 iterations done. 312821144 characters processed.
[BWTIncConstructFromPacked] 90 iterations done. 334797592 characters processed.
[bwt_gen] Finished constructing BWT in 93 iterations.
[bwa_index] 93.39 seconds elapse.
[bwa_index] Update BWT... 0.95 sec
[bwa_index] Pack forward-only FASTA... 0.59 sec
[bwa_index] Construct SA from BWT and Occ... 37.09 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index Dvir.DiNV.combo.fna
[main] Real time: 133.065 sec; CPU: 133.044 sec
```

- Next is to map the trimmed sequences I have to the indexed genome
- The command bwa mem "Align[s] 70bp-1Mbp query sequences with the BWA-MEM algorithm. Briefly, the algorithm works by seeding alignments with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW)"  
  - That basically means aligning a read based off of an exact match in the sequence to the genome, and then extending out off that exact match
- BWA outputs files in SAM format. SAM files are huge, however they are readable format. Most programs take a BAM file, which is the same as a SAM file but much smaller and are not human readable. I will need to convert the SAM output into a BAM file
- To do this, I need a program called [samtools](http://www.htslib.org/doc/samtools.html)
- The command `samtools view` prints all alignment files, and the options -h means to include the header in the output, and -b means that the output will be a BAM file
- The file also needs to be sorted so that it is in order (it will be random otherwise) and this is done with `samtools sort` which will sort the alignments by the leftmost coordinates   
- These commands can all be "piped" together because they use the output of one command as the input of the next. The "|" is the pipe  
`bwa mem Dvir.DiNV.combo.fna KM_1_trim_d.fq.gz KM_2_trim_d.fq.gz | samtools view -hb - | samtools sort - > KM_3.mapped.bam`
- I also wanted to see what the remove PCR duplicates (-rmdup) flag would do at the bam file stage, so I took the bam file I had and used that function on it   
`samtools rmdup KM_3.mapped.bam KM_3_dd.mapped.bam`
- I checked how many lines in the files to see what the removing duplication did   
- Original sequence file:  
`zcat KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz | echo $((`wc -l`/4))`
7150325  
- Trimmed file without deduplication:   
`zcat KM_3_1_trim.fq.gz | echo $((`wc -l`/4))`  
7026662
- Trimmed file with deduplication:  
`zcat KM_3_1_trim_d.fq.gz | echo $((`wc -l`/4))`  
6367634
- Also check the file sizes of the two bam files
- Non-deduplicated file:  
`ls -s KM_3.mapped.bam`  
803976
- Deduplicated file:  
`ls -s KM_3_dd.mapped.bam`  
773408
- I know we looked at the bam files and saw what looked like duplicate sequences, so I want to check those files again to see if they're like that in the deduplicated file. The non-samtools-rmdup bam file was made with fastp-dedup sequences
- View the bam file, looking for DiNV in the 3rd column, not DiNV in the 7th column, and not = in the 7th column. This means one read maps to DiNV, the read pair doesn't map to DiNV, and it also has a pair (= means not paired)  
- Mapped bam:  
`samtools view KM_3.mapped.bam | awk '$3=="DiNV" && $7!="DiNV" && $7!="="'`
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/mapped-bam.png)
- rmdup mapped bam:   
`samtools view KM_3_dd.mapped.bam | awk '$3=="DiNV" && $7!="DiNV" && $7!="="'`
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/rmpdup-mapped-bam.png)
- These produce the same output no matter the file, so it seems like those sequences that look duplicated are not?
