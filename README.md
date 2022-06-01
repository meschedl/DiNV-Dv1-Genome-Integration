# DiNV-Dv1-Genome-Integration

In this repository I investigated whether sequence data from DV-1 cells infected with DiNV showed any evidence of DiNV integrating into the DV-1 genome. 
This seemed possible because the DiNV genome replication rate in DV-1 cells is basically the same as the DV-1 geome replication rate. 

I did some QC on the sequnce data first (only 1 sample forward and reverse) in the [Data_QC](https://github.com/meschedl/DiNV-Dv1-Genome-Integration/tree/main/Data_QC) folder 

Then I mapped the reads to a concatinated DiNV-D. virilis genome, and looked for chimeric reads. These are reads that map both to DiNV and to a virilis chromosome. We were looking for a "half and half" senerio. These reads were investigated and found to just be regions of the reads that map perfectly to both DiNV and virilis. The code for those investegations is in the [Mapping](https://github.com/meschedl/DiNV-Dv1-Genome-Integration/tree/main/Mapping) folder. 

We have concluded that there is no evidence that DiNV integrates into the genome of DV-1 cells during infection. 
