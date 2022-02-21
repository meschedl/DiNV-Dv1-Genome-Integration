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
- Secure copy in the _D. virilis_ genome   
`scp /Users/m741s365/Desktop/chr*  runcklesslab@10.119.46.137:/home/runcklesslab/Maggie`
- Secure copy in the DiNV genome   
`scp /Users/m741s365/Desktop/GCF_004132165.1_ASM413216v1_genomic.fna  runcklesslab@10.119.46.137:/home/runcklesslab/Maggie`
- I want to map my sequences to both the DiNV genome and the _D. virilis_ genome at the same time, so they will have to be concatenated. I also want to give them recognizable header lines. Right now they look like this: `CM017605.2 Drosophila virilis strain 160 chromosome 2, whole genome shotgun sequence` or `NC_040699.1 Drosophila innubila nudivirus isolate DiNV_CH01M, complete genome`
- Only the first word gets read, so I want to change that with nano   
`nano chr2.fna`  
Add in chr2 so that it reads `>chr2 CM017605.2 Drosophila virilis strain 160 chromosome 2, whole genome shotgun sequence`  
Use control (^) O to save and ^X to exit the file
- Also want to change the DiNV header sequence  
`nano GCF_004132165.1_ASM413216v1_genomic.fna`  
Add in DiNV so that it reads `DiNV NC_040699.1 Drosophila innubila nudivirus isolate DiNV_CH01M, complete genome`
- There is also a file unplaced.scaf.fna that I don't know what to do with. It is from the _D. virilis_ genome. When you look at it, it is a bunch of repeats of Ts and Gs, with a few As and no Cs. So for right now I'm not going to use it?

**Step 2: Fastp Sequence Quality Checks**
- [Fastp](https://github.com/OpenGene/fastp) 







 1004 fastp
 1005 fastp -i KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz -I KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_2.fq.gz -o KM_3_1.fq.gz -O KM_3_2.fq.gz
 1006 ls
 1007 head chr2.fna
 1008 nano chr3.fna
 1009 history
