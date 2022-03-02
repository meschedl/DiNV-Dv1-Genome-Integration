## Copying Data and Genomes into Folder and Running Fastp on Data

**Step 1: Access Lab Linux Machine, Create a Directory, and Move my Files**  
- ssh into machine
- Create directory   
`mkdir Maggie`  
`cd Maggie`
- Copy sequence data into directory   
`cp /mnt/data/RawReads/Spring2020Seq/X202SC20021689-Z01-F001/raw_data/KM_3_DNA/KM* .`
- Sequence files are called:
  - KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz
  - KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_2.fq.gz
- Download the DiNV genome and the _D. virilis_ genome to my computer to secure copy to the Linux machine
  - [DiNV genome ASM413216v1](https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=2057187&utm_source=gquery&utm_medium=referral) downloaded as FASTA sequence
  - [_D. virilis_ genome ASM798932v2](https://www.ncbi.nlm.nih.gov/datasets/genomes/?taxon=7244&utm_source=gquery) downloaded as FASTA sequence, we wanted chromosome assembled
- Secure copy in the _D. virilis_ genome. **Use desktop terminal window for this**     
`scp /Users/m741s365/Desktop/chr*  runcklesslab@10.119.46.137:/home/runcklesslab/Maggie`
- Secure copy in the DiNV genome **Use desktop terminal window for this**  
`scp /Users/m741s365/Desktop/GCF_004132165.1_ASM413216v1_genomic.fna  runcklesslab@10.119.46.137:/home/runcklesslab/Maggie`
- I want to map my sequences to both the DiNV genome and the _D. virilis_ genome at the same time, so they will have to be concatenated. I also want to give them recognizable header lines. Right now they look like this: `CM017605.2 Drosophila virilis strain 160 chromosome 2, whole genome shotgun sequence` or `NC_040699.1 Drosophila innubila nudivirus isolate DiNV_CH01M, complete genome`
- Only the first word gets read, so I want to change that with nano   
`nano chr2.fna`  
Add in chr2 so that it reads `>chr2 CM017605.2 Drosophila virilis strain 160 chromosome 2, whole genome shotgun sequence`  
Use control (^) O to save and ^X to exit the file
- Also want to change the DiNV header sequence  
`nano GCF_004132165.1_ASM413216v1_genomic.fna`  
Add in DiNV so that it reads `DiNV NC_040699.1 Drosophila innubila nudivirus isolate DiNV_CH01M, complete genome`
- There is also a file unplaced.scaf.fna that I don't know what to do with. It is from the _D. virilis_ genome. When you look at it, it is a bunch of repeats of Ts and Gs, with a few As and no Cs. It has unplaced as it's first word, so I'll keep that 

**Step 2: Fastp Sequence Quality Checks**
- [Fastp](https://github.com/OpenGene/fastp) is a program that will give you various quality metrics on sequence data, and also allow you to do trimming etc. on the sequences
- The default usage does a couple of things:
  - Removes bases with a lower than 15 phred score
  - Adapter trimming
  - Length filtering where reads shorter than 15bp are removed
  - Gives you an output file that is an html file
- I ran the standard fastp on the data (just one forward and reverse read) to see what it looked like and if any other quality filtering would be needed   
`fastp -i KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz -I KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_2.fq.gz -o KM_3_1.fq.gz -O KM_3_2.fq.gz`
- This outputs the trimmed files as the shorter names KM_3_1.fq.gz and KM_3_2.fq.gz
- To look at the output file, I need to secure copy it to my desktop and open. **Use desktop terminal window for this**   
`scp runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/fastp.html /Users/m741s365/Desktop/Github/DiNV-Dv1-Genome-Integration/Data_QC/`
- [This is what the output file looks like](https://rawcdn.githack.com/meschedl/DiNV-Dv1-Genome-Integration/6b880e89749ed66f0f7a073d77e667b6a61ac702/Data_QC/fastp.html)
- There are a few main things I noticed
  - The quality trails off more for read 2 than read 1 at the end of the sequence
  - The first couple of bases have low scores, so I'll want to trim those, looks like 7 bases for each read. Even though they're still over 30 in score it looks weird to me
  - There is also a lot of noise in the GC content in the first ~15 bases. That should be constant across the read
  - There is a drop in quality score around base 76-77 for both the reads, I'm not sure if there's anything I can do about it. It only drops to a score of 37 which is still really good
- Going to run fastp again and trim the first 15 bases, other than that it probably is good   
`fastp -i KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz -I KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_2.fq.gz --trim_front1 15 --trim_front2 15 --html trim-fastp.html -o KM_3_1_trim.fq.gz -O KM_3_2_trim.fq.gz`
- Here is the output from the program:

```
Read1 before filtering:
total reads: 7150325
total bases: 1072548750
Q20 bases: 1052709214(98.1502%)
Q30 bases: 1030796684(96.1072%)

Read2 before filtering:
total reads: 7150325
total bases: 1072548750
Q20 bases: 1032920148(96.3052%)
Q30 bases: 994062589(92.6823%)

Read1 after filtering:
total reads: 7026662
total bases: 812090153
Q20 bases: 799654639(98.4687%)
Q30 bases: 784993787(96.6634%)

Read2 aftering filtering:
total reads: 7026662
total bases: 812090153
Q20 bases: 787620681(96.9869%)
Q30 bases: 758952060(93.4566%)

Filtering result:
reads passed filter: 14053324
reads failed due to low quality: 247324
reads failed due to too many N: 2
reads failed due to too short: 0
reads with adapter trimmed: 7718326
bases trimmed due to adapters: 273565416

Duplication rate: 11.4592%

Insert size peak (evaluated by paired-end reads): 170
```

- To look at the output file, I need to secure copy it to my desktop and open. **Use desktop terminal window for this**   
`scp runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/trim-fastp.html /Users/m741s365/Desktop/Github/DiNV-Dv1-Genome-Integration/Data_QC/`
- [This is what the output file looks like](https://rawcdn.githack.com/meschedl/DiNV-Dv1-Genome-Integration/6b880e89749ed66f0f7a073d77e667b6a61ac702/Data_QC/trim-fastp.html)
- Move all the QC files and info into a new directory  
`mkdir Seq_QC`  
`mv /home/runcklesslab/Maggie/trim-fastp.html /home/runcklesslab/Maggie/Seq_QC/trim-fastp.html`  
`mv /home/runcklesslab/Maggie/fastp.html /home/runcklesslab/Maggie/Seq_QC/fastp.html`  
`mv /home/runcklesslab/Maggie/fastp.json /home/runcklesslab/Maggie/Seq_QC/fastp.json`  
`mv /home/runcklesslab/Maggie/KM_3_1.fq.gz /home/runcklesslab/Maggie/Seq_QC/KM_3_1.fq.gz`  
`mv /home/runcklesslab/Maggie/KM_3_2.fq.gz /home/runcklesslab/Maggie/Seq_QC/KM_3_2.fq.gz`


**After talking with Rob, we decided that the 11% PCR duplication rate is pretty high, and that I should try to remove PCR duplicates with fastp**
- Again I'll run the program on the original files with the added --dedup deduplication flag (need to make sure version of fastp is past v0.22.0, I had to download this to my folder and run from there)  
`/home/runcklesslab/Maggie/fastp -i KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz -I KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_2.fq.gz --trim_front1 15 --trim_front2 15 --dedup --html trim-dedup-fastp.html -o KM_3_1_trim_d.fq.gz -O KM_3_2_trim_d.fq.gz`
- Copy the file to my computer and send to Github   
`scp runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/trim-dedup-fastp.html /Users/m741s365/Desktop/Github/DiNV-Dv1-Genome-Integration/Data_QC/`
- [Here is the deduplicated and trimmed file](https://rawcdn.githack.com/meschedl/DiNV-Dv1-Genome-Integration/1c01f4b2e49f182a4cf41e7e661d874087462a5d/Data_QC/trim-dedup-fastp.html)
- The problem of the shaky G,C,A, and T content lines doesn't go away.
