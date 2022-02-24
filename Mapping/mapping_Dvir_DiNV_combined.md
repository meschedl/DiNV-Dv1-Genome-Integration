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
`ln -s /home/runcklesslab/Maggie/KM_3_1_trim.fq.gz Mapping/KM_1_trim.fq.gz`  
`ln -s /home/runcklesslab/Maggie/KM_3_2_trim.fq.gz Mapping/KM_2_trim.fq.gz`

**Step 2: Using BWA to Create a Genome Index and Map to the Genome**
- [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) is an alignment tool that is widely used
- First "BWA first needs to construct the FM-index for the reference genome"  
`bwa index Dvir.DiNV.combo.fna`  
- Output from bwa index:   

```
[bwa_index] Pack FASTA... 1.01 sec
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
[bwa_index] 85.44 seconds elapse.
[bwa_index] Update BWT... 0.86 sec
[bwa_index] Pack forward-only FASTA... 0.58 sec
[bwa_index] Construct SA from BWT and Occ... 32.71 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index Dvir.DiNV.combo.fna
[main] Real time: 120.618 sec; CPU: 120.595 sec
```

- Next is to map the trimmed sequences I have to the indexed genome
- The command bwa mem "Align[s] 70bp-1Mbp query sequences with the BWA-MEM algorithm. Briefly, the algorithm works by seeding alignments with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW)"  
  - That basically means aligning a read based off of an exact match in the sequence to the genome, and then extending out off that exact match
- BWA outputs files in SAM format. SAM files are huge, however they are readable format. Most programs take a BAM file, which is the same as a SAM file but much smaller and are not human readable. I will need to convert the SAM output into a BAM file
- To do this, I need a program called [samtools](http://www.htslib.org/doc/samtools.html)
- The command `samtools view` prints all alignment files, and the options -h means to include the header in the output, and -b means that the output will be a BAM file
- The file also needs to be sorted so that it is in order (it will be random otherwise) and this is done with `samtools sort` which will sort the alignments by the leftmost coordinates. And the -o option gives it an output file name   
- These commands can all be "piped" together because they use the output of one command as the input of the next. The "|" is the pipe  
`bwa mem Dvir.DiNV.combo.fna KM_1_trim.fq.gz KM_2_trim.fq.gz | samtools view -h -b | samtools sort -o KM_3.mapped.bam`
