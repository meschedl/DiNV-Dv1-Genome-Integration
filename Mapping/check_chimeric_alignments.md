## Investigating Chimeric Alignments

Assessing the validity of the chimeric alignments to check if they are just regions of the DV-1 and DiNV genome that are overlapping.

For each chimeric alignment, I will go through and pull out the read, and the BAM alignments for the read and align them with the [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Based off of the alignment I will go forward with that chimeric alignment or not.


Starting with **chr3 6251589**  
- In the Linux   
`cd /home/runcklesslab/Maggie/DiNV-DV-1/Mapping`
- Looking for read E00489:560:H7N33CCX2:1:2220:8927:35027
- This read is actually on there twice, probably the forward and reverse read?
- Pull out the read:  
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:2220:8927:35027"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`  
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:2220:8927:35027"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_2.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:2220:8927:35027`
- Used [this website](https://www.bioinformatics.org/sms/rev_comp.html) to make the reverse compliment of the DiNV read
- Copied them to the Clustal page, ended up not using the reverse read
```
>forward
GCGTTCACAATATATGAATTTAAACGTTAAGTTTTTTTATTTCGAATTATTATTTTTTTTATTCAAAAATATATACAATCATTTTCAGATTATTATTGAAGGGAGAGGGGCGACGGATCAGACGGAGCAGGAACTGGCAACAGCCCGGAC
>chr3
GGCGACGGATCAGACGGAGCAGGAACTGGCAACAGCCCGGAC
>DiNV
GAATTTAAACGTTAAGTTTTTTTATTTCGAATTATTATTTTTTTTATTCAAAAATATATACAATCATTTTCAGATTATTATTGAAGGGAGAGGGGCGACGGATCAGACGGAGCAGGAACTGGCAACAGCCCGGAC
```
- Based off of the alignment, they DiNV and Chr3 sections seem to be completely overlapping and this doesn't look like the type of chimeric read we need
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/6251589_align.png)

**chrx 10751602**
- Pull out read E00489:560:H7N33CCX2:1:1103:6431:15619  
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:1103:6431:15619"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:1103:6431:15619`
- This is interesting, this read maps to chr3, chr2, chx, and DiNV. Not looking good!
- Used the chrx alignment and the reverse compliment of the DiNV alignment again, just by looking by eye this seemed like the right thing to do
```
>forward
TTTCATTGGTGCCCAAATAGAAACGTTGACCTTGTCTTGCAGCGATAGGTTGTAAACATATTGCATTTGTATTTAATTTAACTAAACTTTCATATTTTTCACATTGCCATTCCGATGATAATTGTTGTTGTTGCTGTTGTTGTTGTTGTT
>chrx
GATGATAATTGTTGTTGTTGCTGTTGTTGTTGTTGTT
>DiNV
AATAGAAACGTTGACCTTGTCTTGCAGCGATAGGTTGTAAACATATTGCATTTGTATTTAATTTAACTAAACTTTCATATTTTTCACATTGCCATTCCGATGATAATTGTTGTTGTTGCTGTTGTTGTTGTTGTT
```
- Again the two alignments completely overlap, meaning this isn't a half and half read at all
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/10751602-align.png)

**chr3 19526083**
- Pull out read E00489:560:H7N33CCX2:1:1119:30289:8148
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:1119:30289:8148"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:1119:30289:8148`
- This read maps to DiNV, chr3, and chrx...
- I did not do the reverse compliment of either of the sequences because they looked like they would align without it
```
>forward
GATAAATGCCACACGCATGTGGCATGCAACAAAACGATTTTATATGCTTATTATGATTACTGACAGTTGTCGGTTGTTTCTGCTGCTGCTGCTGCTGCTGCTGCTGCGGCTGCTGCTGCTGTGCTCTGCTGTGGTGGCTGTCTCTTATAC
>chr3
GATAAATGCCACACGCATGTGGCATGCAACAAAACGATTTTATATGCTTATTATGATTACTGACAGTTGTCGGTTGTTTCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTACTGATGT
>DiNV
TGCTGCTGCTGCTGCTGCTGCTGCTACTGATGT
```
- And that is correct, the alignments from DiNV and Chr3 are completely overlapping
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/19526083_align.png)

**chrx 17669245**
- Pull out read E00489:560:H7N33CCX2:1:1205:6969:4983
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:1205:6969:4983"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:1205:6969:4983`
- Used the reverse compliment of the DiNV sequence
```
>forward
GCCATTCCGATGATAATTGTTGTTGTTGCTGTTGTTGTTGTTGTTGTTGTACTTCTTTTTTTGATTCTATTAAAGACGACGATGACGATGATGATACTGATGTAGAAAATAATGAAGAATTCAAAACATCTAATAAACTTCCATTAGCAT
>chrx
ATTGTTGTTGTTGCTGTTGTTGTTGTTGTTGTTGTAC
>DiNV
GTTGTTGCTGTTGTTGTTGTTGTTGTTGTACTTCTTTTTTTGATTCTATTAAAGACGACGATGACGATGATGATACTGATGTAGAAAATAATGAAGAATTCAAAACATCTAATAAACTTCCATTAGCATTTTTAA
```
- The alignments from chrx and DiNV are completely overlapping
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/17669245_align.png)

**chr2 32783122**
- Pull out read E00489:560:H7N33CCX2:1:1103:6431:15619
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:1106:12063:16164"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:1106:12063:16164`
- This read maps to both chr2 and chrx, as well as DiNV
```
>forward
CTTTAATAGAATCAAAAAAAGAAGTACAACAACAACAACAACAACAGCAACAACAACAATTATCATCGGAATGGCAATGTGAAAAATATGAAAGTTTAGTTAAATTAAATACAAATGCAATATGTTTACAACCTATCGCTGCAAGACAAG
>chr2
TACAACAACAACAACAACAACAGCAACAACAACAA
>DiNV
AGAATCAAAAAAAGAAGTACAACAACAACAACAACAACAGCAACAACAACAATTATCATCGGAATGGCAATGTGAAAAATATGAAAGTTTAGTTAAATTAAATACAAATGCAATATGTTTACAACCTATCGCTGC
```
- The alignments from chr2 and DiNV are completely overlapping
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/32783122_align.png)

**DiNV 39583**
- Pull out read E00489:560:H7N33CCX2:1:1119:30289:8148
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:1119:30289:8148"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:1119:30289:8148`
- This read maps to chr3, chrx, and DiNV
```
>forward
GATAAATGCCACACGCATGTGGCATGCAACAAAACGATTTTATATGCTTATTATGATTACTGACAGTTGTCGGTTGTTTCTGCTGCTGCTGCTGCTGCTGCTGCTGCGGCTGCTGCTGCTGTGCTCTGCTGTGGTGGCTGTCTCTTATAC
>chr3
GATAAATGCCACACGCATGTGGCATGCAACAAAACGATTTTATATGCTTATTATGATTACTGACAGTTGTCGGTTGTTTCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTACTGATGT
>DiNV
TGCTGCTGCTGCTGCTGCTGCTGCTACTGATGT
```
- The alignments from chr3 and DiNV are completely overlapping
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/39583_align.png)

**DiNV 52112**
- Pull out read E00489:560:H7N33CCX2:1:2220:8927:35027
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:2220:8927:35027"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:2220:8927:35027`
- I used the reverse compliment of the DiNV sequence
```
>forward
GCGTTCACAATATATGAATTTAAACGTTAAGTTTTTTTATTTCGAATTATTATTTTTTTTATTCAAAAATATATACAATCATTTTCAGATTATTATTGAAGGGAGAGGGGCGACGGATCAGACGGAGCAGGAACTGGCAACAGCCCGGAC
>chr3
GGCGACGGATCAGACGGAGCAGGAACTGGCAACAGCCCGGACGAACCTGC
>DiNV
GAATTTAAACGTTAAGTTTTTTTATTTCGAATTATTATTTTTTTTATTCAAAAATATATACAATCATTTTCAGATTATTATTGAAGGGAGAGGGGCGACGGATCAGACGGAGCAGGAACTGGCAACAGCCCGGAC
```
- The alignments are completely overlapping
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/52112_align.png)

**DiNV 132752 (there are actually multiple chimeric alignments in this region so I'm going to do a couple of reads)**
- Pull out read E00489:560:H7N33CCX2:1:1106:12063:16164
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:1106:12063:16164"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:1106:12063:16164`
- This read maps to chrx, chr2, and DiNV
```
>forward
CTTTAATAGAATCAAAAAAAGAAGTACAACAACAACAACAACAACAGCAACAACAACAATTATCATCGGAATGGCAATGTGAAAAATATGAAAGTTTAGTTAAATTAAATACAAATGCAATATGTTTACAACCTATCGCTGCAAGACAAG
>chr2
TACAACAACAACAACAACAACAGCAACAACAACAAT
>DiNV
AGAATCAAAAAAAGAAGTACAACAACAACAACAACAACAGCAACAACAACAATTATCATCGGAATGGCAATGTGAAAAATATGAAAGTTTAGTTAAATTAAATACAAATGCAATATGTTTACAACCTATCGCTGC
```
- The alignments for DiNV and chr2 are completely overlapping
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/132752_1_align.png)
- Pull out read E00489:560:H7N33CCX2:1:2202:28686:65775
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:2202:28686:65775"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:2202:28686:65775`
- Used reverse compliment of the DiNV alignment
```
>forward
GTTTCATTGGTGCCCAAATAGAAACGTTGACCTTGTCTTGCAGCGATAGGTTGTAAACATATTGCATTTGTATTTAATTTAACTAAACTTTCATATTTTTCACATTGCCATTCCGATGATAATTGTTGTTGTTGCTGTTGTTGTTGTTGT
>chrx
GATGATAATTGTTGTTGTTGCTGTTGTTGTTGTTGT
>DiNV
AAATAGAAACGTTGACCTTGTCTTGCAGCGATAGGTTGTAAACATATTGCATTTGTATTTAATTTAACTAAACTTTCATATTTTTCACATTGCCATTCCGATGATAATTGTTGTTGTTGCTGTTGTTGTTGTTGT
```
- The alignments for DiNV and chrx completely overlap
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/132752_2_align.png)
- Pull out read E00489:560:H7N33CCX2:1:1205:6969:4983
`bioawk -c fastx '$name=="E00489:560:H7N33CCX2:1:1205:6969:4983"' ../KM_3_DNA_CKDL200150169-1a-N705-N505_H7N33CCX2_L1_1.fq.gz`
- Pull out the alignments:  
`samtools view KM_3.mappedd_F4.bam | grep E00489:560:H7N33CCX2:1:1205:6969:4983`
- Used the reverse compliment of the DiNV sequence
```
>forward
GCCATTCCGATGATAATTGTTGTTGTTGCTGTTGTTGTTGTTGTTGTTGTACTTCTTTTTTTGATTCTATTAAAGACGACGATGACGATGATGATACTGATGTAGAAAATAATGAAGAATTCAAAACATCTAATAAACTTCCATTAGCAT
>chrx
ATTGTTGTTGTTGCTGTTGTTGTTGTTGTTGTTGTAC
>DiNV
ATTGTTGTTGTTGCTGTTGTTGTTGTTGTTGTTGTACTTCTTTTTTTGATTCTATTAAAGACGACGATGACGATGATGATACTGATGTAGAAAATAATGAAGAATTCAAAACATCTAATAAACTTCCATTAGCAT
```
- The chrx and DiNV sequences are completely overlapping
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/132752_3_align.png)


Based off of these, there are no "good" chimeric alignments where half the read is DV-1 and the other half of the read is DiNV and are completely separate. That was what we are looking for. There are more chimeric alignments that I didn't look at, either I ruled out because they seemed very repetitive (and these seemed this way too) or they only had 1 read coverage. We are not able to find any evidence that DiNV integrates into DV-1 cell genome with this method.
