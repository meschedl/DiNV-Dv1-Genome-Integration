## Generating Primers to Try Amplifying Chimeric Junctions

My goal is to generate primers that span some of the chimeric junctions in the sequencing of DV-1 cells infected with DiNV. I will need to generate primers that amplify the Dv-1 region around the chimeric alignment and the DiNV region around the chimeric alignment. Then, a combination of those two would in theory amplify in the integrated sample (see image).

![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/primer-image.png)

**Starting with chrx position 10751602**
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/10751602-all.png)
- There are 5 reads here that have that chimeric alignment
- I'm going to pick the 1st read E00489:560:H7N33CCX2:1:1103:6431:15619 and see what is going on in the bam file

`cd /home/runcklesslab/Maggie/DiNV-DV-1/Mapping/chimeric_reads`  
`grep E00489:560:H7N33CCX2:1:1103:6431:15619 all_chimeric_KM3.sam`
- Looking at this, the read maps to DiNV, chr X, and Chr 2
- I only care about DiNV and Chr X here, but it looks like the bases that map to chr X are very few: GATGATAATTGTTGTTGTTGCTGTTGTTGTTGTTGTT
- The reality of generating a primer there seems impossible...
- I will try anyways
- This is at position 10751602 in DV-1 genome on the X chromosome
- If I make a reverse primer from this sequence, it's:
  - ACAACAACAACAACAACAACAACA
  - GC% 33.3 (too low)
  - Tm 59.1
- That primer sequence seems so repetitive to me I can't imagine it having specific binding at all. I'm going to try to look for a different chimeric alignment to try this on

**chr3 6251589**
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/6251589-all.png)
- Look at read E00489:560:H7N33CCX2:1:2220:8927:35027
- Looks like it maps to ~50 bases in chr3 and 55M1D22M1D15M to DiNV
- Things are clipped on the left so I'm going to make a reverse primer for the chimeric alignment
- Reverse ch3 primer attempt 1:
  - CTGTTGCCAGTTCCTGCTCC
  - GC% 60.0
  - Tm 61.2
- Reverse ch3 primer attempt 2:
  - GTTCCTGCTCCGTCTGATCC
  - GC% 60.0
  - Tm 60.1
- These may have too many Cs but that is just the nature of the sequence, there isn't much room to create a primer in there?
- Forward primers for chr3 region
- Forward primer 1
  - CAAGTCAAAATAGAGGCTACCAGA
  - GC% 41.7
  - Tm 58.3
  - Starts at 6,251,307
- Forward primer 2
  - CCAATTGATGACTCGGTTCGG
  - GC% 52.4
  - Tm 59.3
  - Starts at 6,251,284

![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/chr3-6251284-primers.png)
- Now to find the DiNV region
- This is at 52112 in DiNV
