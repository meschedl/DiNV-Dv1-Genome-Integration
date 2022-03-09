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
- The command `samtools view` prints all alignment files, and the options -h means to include the header in the output
- The file also needs to be sorted so that it is in order (it will be random otherwise) and this is done with `samtools sort` which will sort the alignments by the leftmost coordinates   
- These commands can all be "piped" together because they use the output of one command as the input of the next. The "|" is the pipe  
`bwa mem Dvir.DiNV.combo.fna KM_1_trim_d.fq.gz KM_2_trim_d.fq.gz | samtools view -h | samtools sort - > KM_3.mapped-h.bam`
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
- These produce the same output no matter the file, so it seems like those sequences that look duplicated are not

**Filtering BAM File and Looking for Chimeric Reads**
- From what I can tell with searching the internet, is that the MAPQ column (which is supposed to be an indication of mapping quality) does not really tell you anything when using BWA-mem
- One thing I can do easily is remove reads that aren't mapped. I can do this by using the -F flag in samtools view, which will not output alignments with the specified FLAG. If I do -F 4 that will not output alignments where the sequence is not mapped. I am going to do this with the de-duplicated file   
`samtools view -h -b -F 4 KM_3_dd.mapped.bam > KM_3.mappedd_F4.bam`
- Sanity check that the files are different  
`ls -s KM_3.mapped_dd.bam`  
773408  
`ls -s KM_3.mappedd_F4.bam`  
762144
- From this point, I am not sure what else I can do to filter the BAM file to be only alignments that we trust. I am going to move on to looking for the indicators of DiNV integration: either reads that map to a chromosome and the pair maps to DiNV, or chimeric reads that map to both a chromosome and DiNV
- Chimeric reads are indicated in the SAM format in the 12th column with an SA followed by other words. The following text is rname, position, strand, CIGAR string, mapQ, NM (number of differences between the sequence and the reference). I am pretty sure if the read is not chimeric there is not an SA flag in the 12th column. So what I need to do is search the BAM file for only lines that include an SA  
`samtools view -h KM_3.mappedd_F4.bam | grep 'SA:' > KM_3_SA.bam`
- Sanity check the file sizes between these files  
`ls -s KM_3.mappedd_F4.bam`   
762144  
`ls -s KM_3_SA.bam`  
29448
- Is there a difference between the files that grep SA: and SA:Z:?  
 `samtools view -h KM_3.mappedd_F4.bam | grep 'SA:Z:' > KM_3_SAZ.bam`  
- Check size:
 `ls -s KM_3_SAZ.bam`  
 29448
- Nope, same size. This makes further grep-ing/separating out easier I think. Why? Because an example of the SA flag looks like this `SA:Z:chr6,16979,-,41M94S,0,0;` where it has the reference chromosome, the position, the strand, the CIGAR string, the MAPQ and the NM. I will want to separate out who has a chimeric read to specific chromosomes/DiNV
- See if this works, first I will separate out this file to everything that has the initial mapping to chr2 , then pipe it into grep where i will search for all the lines where the chimeric read also maps to DiNV and save that into a new file
`samtools view KM_3_SAZ.bam | awk '$3=="chr2"' | grep 'SA:Z:DiNV' > chim_chr2_DiNV.bam`
- this doesn't work because samtools view wants the header lines in the bam file! why!! Ok what I can do is copy and paste the header lines into the bam file, I loose them when I grep because it doesn't use those lines. Maybe I could also look for those lines too??
`samtools view -h KM_3.mappedd_F4.bam | less`
- These are the header lines:
```
@HD     VN:1.6  SO:coordinate
@SQ     SN:chr2 LN:38438298
@SQ     SN:chr3 LN:27616680
@SQ     SN:chr4 LN:31075311
@SQ     SN:chr5 LN:27902728
@SQ     SN:chr6 LN:2270151
@SQ     SN:chrx LN:38193566
@SQ     SN:DiNV LN:155555
@SQ     SN:unplaced     LN:16908
@SQ     SN:VNHH02000019.1       LN:560307
@SQ     SN:VNHH02000043.1       LN:331422
@SQ     SN:VNHH02000047.1       LN:234219
@SQ     SN:VNHH02000048.1       LN:20672
@SQ     SN:VNHH02000050.1       LN:275122
@SQ     SN:VNHH02000051.1       LN:29855
@SQ     SN:VNHH02000054.1       LN:73214
@SQ     SN:VNHH02000055.1       LN:19572
@SQ     SN:VNHH02000056.1       LN:52186
@SQ     SN:VNHH02000082.1       LN:47289
@SQ     SN:VNHH02000089.1       LN:52137
@SQ     SN:VNHH02000090.1       LN:34989
@SQ     SN:VNHH02000095.1       LN:19220
@SQ     SN:VNHH02000097.1       LN:29477
@SQ     SN:VNHH02000102.1       LN:56772
@SQ     SN:VNHH02000105.1       LN:30176
@SQ     SN:VNHH02000107.1       LN:21993
@SQ     SN:VNHH02000108.1       LN:19534
@SQ     SN:VNHH02000111.1       LN:26096
@SQ     SN:VNHH02000112.1       LN:17191
@SQ     SN:VNHH02000113.1       LN:49501
@SQ     SN:VNHH02000116.1       LN:21351
@SQ     SN:VNHH02000127.1       LN:25842
@SQ     SN:VNHH02000150.1       LN:521567
@SQ     SN:VNHH02000151.1       LN:86477
@SQ     SN:VNHH02000171.1       LN:104324
@SQ     SN:VNHH02000177.1       LN:231566
@SQ     SN:VNHH02000181.1       LN:21731
@SQ     SN:VNHH02000187.1       LN:442058
@SQ     SN:VNHH02000188.1       LN:321154
@SQ     SN:VNHH02000191.1       LN:107803
@SQ     SN:VNHH02000195.1       LN:27877
@SQ     SN:VNHH02000196.1       LN:18248
@SQ     SN:VNHH02000197.1       LN:24842
@SQ     SN:VNHH02000198.1       LN:144193
@SQ     SN:VNHH02000201.1       LN:110897
@SQ     SN:VNHH02000207.1       LN:17506
@SQ     SN:VNHH02000208.1       LN:29723
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem Dvir.DiNV.combo.fna KM_1_trim_d.fq.gz KM_2_trim_d.fq.gz
@PG     ID:samtools     PN:samtools     PP:bwa  VN:1.13 CL:samtools view -hb -
@PG     ID:samtools.1   PN:samtools     PP:samtools     VN:1.13 CL:samtools sort -
@PG     ID:samtools.2   PN:samtools     PP:samtools.1   VN:1.13 CL:samtools view -h -b -F 4 KM_3_dd.mapped.bam
@PG     ID:samtools.3   PN:samtools     PP:samtools.2   VN:1.13 CL:samtools view -h KM_3.mappedd_F4.bam
```
- All of those VNHH strings are contigs/scaffolds that are in the unplaced.saffolds.fna file from downloading the genome. I am still not sure if I am going to use those
- what if I tried keeping the header lines in when searching for SA
samtools view -h KM_3.mappedd_F4.bam | awk /"@HD"|"@SQ"|"@PG"|"SA:Z"/ > KM_3_SAZh.bam

made file that'c called grep.text that was just a list:
```
@HD
@SQ
@PG
SA:Z
```


samtools view -h KM_3.mappedd_F4.bam | grep -f grep.txt > KM_3_SAZh.bam

this worked to keep the header lines in the grepd bam file 
