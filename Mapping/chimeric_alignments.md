# Looking for Chimeric Reads = Reads that Map Both to DiNV and D. virilis

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
- What if I tried keeping the header lines in when searching for SA, I can do this my making grep look for a list of things
- Made file that's called grep.text that was just a list:
```
@HD
@SQ
@PG
SA:Z
```
- And then remade the SA:Z file to include the header:  
`samtools view -h KM_3.mappedd_F4.bam | grep -f grep.txt > KM_3_SAZh.bam`
- This worked to keep the header lines in the grepd bam file
- Actually, this is no longer a bam file, because with grep it makes it readable, so I no longer have the issue of needing a header for samtools view to look at it. So now I can subset this file to be just ones that have DiNV as the main alignment, and I no longer need to use samtools view    
`awk '$3=="DiNV"' KM_3_SAZh.bam > KM_3_SAZh_DiNV.bam`
- Now to separate this out to the lines that have SA:Z:chr something. Start with chr2 to get all the lines that have a chimeric alingment to chr2 with DiNV as the main aligment   
`grep 'SA:Z:chr2' KM_3_SAZh_DiNV.bam > KM_3_SAZh_DiNV_chr2.bam`
- Repeat for all the other chromosomes (not including the unplaced scaffold ones)  
`grep 'SA:Z:chr3' KM_3_SAZh_DiNV.bam > KM_3_SAZh_DiNV_chr3.bam`  
`grep 'SA:Z:chr4' KM_3_SAZh_DiNV.bam > KM_3_SAZh_DiNV_chr4.bam`   
`grep 'SA:Z:chr5' KM_3_SAZh_DiNV.bam > KM_3_SAZh_DiNV_chr5.bam`   
`grep 'SA:Z:chr6' KM_3_SAZh_DiNV.bam > KM_3_SAZh_DiNV_chr6.bam`  
`grep 'SA:Z:chrx' KM_3_SAZh_DiNV.bam > KM_3_SAZh_DiNV_chrx.bam`
- Now I want to do the opposite, where I search for all the chromosomes as the main alignment and DiNV as the chimeric one  
`grep 'SA:Z:DiNV' KM_3_SAZh.bam > KM_3_SAZh_chimDiNV.bam`
- Then this file can be separated out into the various chromosomes as the main alignments  
`awk '$3=="chr2"' KM_3_SAZh_chimDiNV.bam > KM_3_SAZh_chr2_chimDiNV.bam`  
`awk '$3=="chr3"' KM_3_SAZh_chimDiNV.bam > KM_3_SAZh_chr3_chimDiNV.bam`  
`awk '$3=="chr4"' KM_3_SAZh_chimDiNV.bam > KM_3_SAZh_chr4_chimDiNV.bam`  
`awk '$3=="chr5"' KM_3_SAZh_chimDiNV.bam > KM_3_SAZh_chr5_chimDiNV.bam`   
`awk '$3=="chr6"' KM_3_SAZh_chimDiNV.bam > KM_3_SAZh_chr6_chimDiNV.bam`  
`awk '$3=="chrx"' KM_3_SAZh_chimDiNV.bam > KM_3_SAZh_chrx_chimDiNV.bam`  
- Made a new directory called chimeric_reads and put all of the KM_3_SAZh files in there
- Talking to Rob, he does want all the unplaced contigs/scaffolds in the analysis, so I'm going to have to go back and separate out each one of those as well
- Separate out the lines that have VNHH as the chimeric alignment from the file that has DiNV as the main alignment. I'm going to try to make this all into one file instead of doing all the separate VNHHs:  
`grep 'SA:Z:VNHH*' KM_3_SAZh_DiNV.bam > KM_3_SAZh_DiNV_VNHH.bam`
- This worked, and there is only 1 of these!
- Now to search for lines that have VNHH* as the main alignment and DiNV as the chimeric one:   
`awk '$3=="VNHH*"' KM_3_SAZh_chimDiNV.bam > KM_3_SAZh_VNHH_chimDiNV.bam`
- There are none that are in this one
- Also want to do the same thing for the "unplaced" scaffold/contig   
`grep 'SA:Z:unplaced' KM_3_SAZh_DiNV.bam > KM_3_SAZh_DiNV_unplaced.bam`  
`awk '$3=="unplaced"' KM_3_SAZh_chimDiNV.bam > KM_3_SAZh_unplaced_chimDiNV.bam`
- Neither of those have anything in them. I am pretty sure I ran the code right for these, and it just really is correct that there aren't barely any chimeric reads for the unplaced sequences
- Now I want to look at these alignments a little bit, do they look really repetitive or do they look like "diverse" sequences
- Look at KM_3_SAZh_DiNV_chr2.bam  
`less KM_3_SAZh_DiNV_chr2.bam`
- There are 3 alignments in this file, their sequences are:
```
> CTGATGGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGTAATTGTGCTTGTTGTACTGGTAATTGCACAACGGGTGTGGATGTGCTAATGGCACCCAATGATTGCATATTAGTGTTTCCATTAGAGGGTTTTACAT  
> AGAATCAAAAAAAGAAGTACAACAACAACAACAACAACAGCAACAACAACAATTATCATCGGAATGGCAATGTGAAAAATATGAAAGTTTAGTTAAATTAAATACAAATGCAATATGTTTACAACCTATCGCTGC  
> AAAAAAAGAAGTACAACAACAACAACAACAACAGCAACAACAACAATTATCATCGGAATGGCAATGTGAAAAATATGAAAGTTTAGTTAAATTAAATACAAATGCAATATGTTTACAACCTATCGCTGCAAGACA
```
- These look somewhat repetitive, but I find it hard to tell if they look like junk or not
- I'm going to look at the CIGAR string, both for the main alignment and for the chimeric alignment. For the first sequence:
  - 24S111M
  - 3S35M97S
  - S means soft clip and M means match
- This means that 111 bases match to DiNV, and 35 bases match to chr2
- Without knowing what those 35 bases are that match to chr2, it's hard to say if it's a realistic chimeric read
- Another question to ask is if the reciprocal file KM_3_SAZh_chr2_chimDiNV.bam has the same number of alignments   
`less KM_3_SAZh_chr2_chimDiNV.bam`
- It does not, there are 6 alignments in here which is interesting. I would have thought that there would be 3
- These are all smaller as well, and look pretty repetitive
```
> CAACAACAACGACAACAGCAACAACAACAA
> TAAAGGGTGCAAATAGAAAAATCTGAATAAAGGGAG
> ATGGTTGTTGTTGTTGTTGTTGTTGTTGTTGTTGT
> TACAACAACAACAACAACAACAGCAACAACAACAAT
> TACAACAACAACAACAACAACAGCAACAACAACAAT
> TACAACAACAACAACAACAACAGCAACAACAACAAT
```
- Those last three are all the same, likely they are just repetitive areas in chr2 that the same sequence mapped to, they have different positions. Although I'm pretty sure they are all different reads because the QNAME is different
- My guess is that these are not "reliable" chimeric reads, mostly based off of the repetitive judgement. Unfortunately that's not a number or something concrete that I can be definitive about
- My best guess on what to do is make a spreadsheet of all the chimeric reads and investigate if they are repetitive, and what the matches are. And I will probably have to go through this by eye, but that's ok because there aren't very many chimeric reads
- I went through every chimeric read by hand and classified them as maybe chimeric or probably not chimeric (based on a really repetitive sequence or a very small amount of the read mapping as chimeric - those I said weren't chimeric) [here](https://github.com/meschedl/DiNV-Dv1-Genome-Integration/blob/main/Mapping/chimeric_alignments_stats.csv)
- Then I made this spreadsheet into a text file in the Linux by using `nano chimeric_stats.txt`
- Then I want to separate out all the lines that say maybe chimeric  
`awk '$8=="maybe"' chimeric_stats.txt > chimeric_stats_maybe.txt`
- Now I want to separate out the first line of all of these, which is the indicator of the read  
`awk '{print$1}' chimeric_stats_maybe.txt > maybe_chimeric_readnames.txt`
- Now I can use this list of read names to separate out my original chimeric reads "BAM" file into just these reads - although it's not a BAM file anymore. If I add the header back in it could read as a SAM file? First I'll add in the header lines to the search document  
`nano maybe_chimeric_readnames.txt` so it has `@HD @SQ @PG` as search parameters  
- Then search the "bam" file for those lines  
`grep -f maybe_chimeric_readnames.txt KM_3_SAZh.bam > maybe_chimeric_KM3.sam`
- This file might be able to be read by IGV? But I can actually put it into a BAM file to to make sure, and IGV prefers a BAM file format  
`samtools sort maybe_chimeric_KM3.sam -o maybe_chimeric_KM3.bam`
- This is a real bam file
- IGV also needs an associated index file with this BAM   
`samtools index -b maybe_chimeric_KM3.bam`
- This creates a .bai file that IGV needs
- Now I want to copy these to my desktop to use IGV on
- Copy combined reference genome to my desktop  
`scp runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/Dvir.DiNV.combo.fa /Users/m741s365/Desktop/`  
`scp runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/Mapping/chimeric_reads/maybe_chimeric_KM3.bam /Users/m741s365/Desktop`  
`scp runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/Mapping/chimeric_reads/maybe_chimeric_KM3.bam.bai /Users/m741s365/Desktop/`
- View in IGV - this is really hard to do because the actual alignments are really small compared to the concatenated genome, and there are so few of them! I tried at first by searching a position in DiNV on the list: 132752 which got me to a region with a lot of overlap in reads that say they are chimeric. This might be a good sign? Looking at the BAM file, there are 10 "alignments" or "reads" that have this position, 8 of them have a chimeric alignment to chrx and 2 of them to chr2. So, maybe this is a sign that it's not true chimeric read because it maps to multiple places in the _D. virilis_ genome? Or it could have inserted into 2 places?
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/132752.png)
- What if I search in IGV all the positions of chimeric reads and image each of those? Want to print both the 3rd column (where it aligned to) and the 4th column (the position)   
`samtools view maybe_chimeric_KM3.bam | awk '{print$3,$4}'`
```
chr2 32783122
chr2 32783122
chr3 6251589
chr3 6251589
chr3 19526083
chr3 19526083
chr3 19526083
chr3 19526097
chr3 19526110
chr3 19526112
chr4 16546625
chr5 1214254
chr5 2877488
chr5 11902161
chrx 4450716
chrx 6782315
chrx 10751602
chrx 10751602
chrx 10751602
chrx 10751602
chrx 10751602
chrx 10751602
chrx 14058146
chrx 14532465
chrx 17669245
chrx 17669246
chrx 19262801
chrx 26750707
chrx 26750707
chrx 27176931
chrx 27639332
chrx 27759345
chrx 32761825
DiNV 1407
DiNV 39583
DiNV 39583
DiNV 39583
DiNV 39583
DiNV 39583
DiNV 44028
DiNV 49076
DiNV 52112
DiNV 52112
DiNV 69463
DiNV 73568
DiNV 84141
DiNV 91291
DiNV 104762
DiNV 109107
DiNV 132635
DiNV 132637
DiNV 132641
DiNV 132703
DiNV 132712
DiNV 132752
DiNV 132752
DiNV 132752
DiNV 132752
DiNV 132752
DiNV 132752
DiNV 132752
DiNV 132752
DiNV 132752
DiNV 132752
DiNV 152624
VNHH02000150.1 78439
```
- Now to look at what all of these look like in IGV

**chr2 32783122**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/32783122.png)
**chr3 6251589**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/6251589.png)
**chr3 19526083**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/19526083.png)
**chr4 16546625**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/16546625.png)
**chr5 1214254**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/1214254.png)
**chr5 2877488**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/2877488.png)
**chr5 11902161**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/11902161.png)
**chrx 4450716**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/4450716.png)
**chrx 6782315**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/6782315.png)
**chrx 10751602**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/10751602.png)
**chrx 14058146**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/14058146.png)
**chrx 14532465**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/14532465.png)
**chrx 17669245**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/17669245.png)
**chrx 19262801**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/19262801.png)
**chrx 26750707**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/26750707.png)
**chrx 27176931**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/27176931.png)
**chrx 27639332**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/27639332.png)
**chrx 27759345**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/27759345.png)
**chrx 32761825**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/32761825.png)
**DiNV 1407**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/1407.png)
**DiNV 39583**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/39583.png)
**DiNV 44028**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/44028.png)
**DiNV 49076**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/49076.png)
**DiNV 52112**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/52112.png)
**DiNV 69463**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/69463.png)
**DiNV 73568**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/73568.png)
**DiNV 84141**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/84141.png)
**DiNV 91291**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/91291.png)
**DiNV 104762**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/104762.png)
**DiNV 109107**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/109107.png)
**DiNV 152624**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/152624.png)
**VNHH02000150.1 78439**  
![](https://raw.githubusercontent.com/meschedl/DiNV-Dv1-Genome-Integration/main/images/IGV_photos/78439.png)

- Problem is that I'm not able to see where the chimeric read pair maps, the BAM file gives me the position of the chimeric alignment, but not the sequence. So maybe the shorter reads that I thought were not "likely" chimeric reads could be the matching ones. I need to match the main alignment position to the SA positions on other reads (hopefully)
- Am I able to look at all chimeric reads and see if the "pairs" are in that file?  
`awk '{print$1}' chimeric_stats.txt > all_chimeric_readnames.txt`  
`grep -f all_chimeric_readnames.txt KM_3_SAZh.bam > all_chimeric_KM3.sam`
`awk '{print$3,$4,$7,$17}' all_chimeric_KM3.sam > all_chimeric_positions.sam`
- This works, but the SA field is pretty messy, I'm going to have to separate it out into tabs: `SA:Z:DiNV,132752,+,47S88M,60,0;`  
`sed 's/,/\t/g' all_chimeric_positions.sam > all_chimeric_positions_t.sam`  
`sed 's/:/\t/g' all_chimeric_positions_t.sam > all_chimeric_positions_tt.sam`  
`awk '{print$1,$2,$3,$6,$7}' all_chimeric_positions_tt.sam > all_chimeric_positions_positions.sam`  
- I want to add column names to this because I might look at this in R to try to figure out if there are any matches of chimerics  
`nano all_chimeric_positions_positions.sam`  
Add `main_alingment main_position mapped_pair  chim_alignment  chim_position`
- Copy this to my desktop to work in R  
`scp runcklesslab@10.119.46.137:/home/runcklesslab/Maggie/Mapping/chimeric_reads/all_chimeric_positions_positions.sam /Users/m741s365/Desktop/`
- This didn't open properly in excel when I saved it as a csv but I can separate out the columns in excel
- Brought this into R to try to see which alignments had a matching position for the chimeric read, looks like they all do
- I think it would be better to have these with all of the read/alignment names as well so I'm going to remake the all_chimeric_positions.sam file to include the 1st column   
`awk '{print$1,$3,$4,$7,$17}' all_chimeric_KM3.sam > all_chimeric_positions.sam`
- Hmm now this is interesting because the read names look like this: `E00489:560:H7N33CCX2:1:1123:12824:3770`, which has a bunch of : in it, which I want to separate out later to get a clear column for the chromosome name of the chimeric read, so actually I don't want to keep these in here, actually just separate them out and paste them into the csv file when I save it to my desktop. There is probably a way to do this in bash but it's too easy to do it this way
- So I am going to go back to the other `all_chimeric_positions.sam`  
`awk '{print$3,$4,$7,$17}' all_chimeric_KM3.sam > all_chimeric_positions.sam`
- So that if I re-do code it'll work right
- And separate out into a separate file just the read "names"  
`awk '{print$1}' all_chimeric_KM3.sam > chimeric_names.txt`
- And copy and paste these into the csv file on my desktop as a new row


I want to bring maybe_chimeric_KM3.sam into R and see if all the positions have a match. The readnames between matches should all be the same, because it's the same read that is chimeric! So trying to match up the read names in a separate column doesn't mean anything .
- I don't need the header lines in maybe_chimeric_KM3.sam so I want to get rid of them
`grep 'E00*' maybe_chimeric_KM3.sam > maybe_chimeric_KM3.txt`
- Now do all the separating out for R
`awk '{print$3,$4,$7,$17}' maybe_chimeric_KM3.txt | sed 's/,/\t/g' | sed 's/:/\t/g' | awk '{print$1,$2,$3,$6,$7}' >  maybe_chimeric_positions.csv`
