# Align reads to the Genome

There are several commonly used align/mapping programs for aligning RNAseq read to the genome.
Hisat2 is a successor of Tophat2 is suitable for align RNAseq reads to the genome.
> - ### Align reads to a reference genome using HISAT2
### [Hisat2 manual](http://daehwankimlab.github.io/hisat2/manual/)

1. Load hisat 2.2.1 module
```  
module load hisat2/2.2.1
```  
2. Creating HISAT2 indexes

First, the genome needed to be indexed before using Hisat2-build or created using the index generatering script:make_mm10.sh
```
cd /scratch/bbash/jl19/Data_QC/

# create a diretory /mm10

mkdir mm10

# Go to the mm10/ direcotry
cd /scratch/bbash/jl19/Data_QC/mm10/

# copy the genome.fa in the mm10 directory

cp /data/bioinf/Teaching/2022_NGS_Course/Data_QC/RNA-Seq-GSE116583/mm10/genome.fa ./

# Creating index

hisat2-build -p 16 genome.fa genome
```
```
[jl19@spectre14 mm10]$ hisat2-build -p 16 genome.fa genome
Settings:
  Output files: "genome.*.ht2"
  Line rate: 6 (line is 64 bytes)
  Lines per side: 1 (side is 64 bytes)
  Offset rate: 4 (one in 16)
  FTable chars: 10
  Strings: unpacked
  Local offset rate: 3 (one in 8)
  Local fTable chars: 6
  Local sequence length: 57344
  Local sequence overlap between two consecutive indexes: 1024
  Endianness: little
  Actual local endianness: little
  Sanity checking: disabled
  Assertions: disabled
  Random seed: 0
  Sizeofs: void*:8, int:4, long:8, size_t:8
Input files DNA, FASTA:
  genome.fa
Reading reference sizes
  Time reading reference sizes: 00:00:21
Calculating joined length
Writing header
Reserving space for joined string
Joining reference sequences
  Time to join reference sequences: 00:00:17
  Time to read SNPs and splice sites: 00:00:00
Using parameters --bmax 31087307 --dcv 1024
  Doing ahead-of-time memory usage test
  Passed!  Constructing with these parameters: --bmax 31087307 --dcv 1024
Constructing suffix-array element generator
Building DifferenceCoverSample
  Building sPrime
  Building sPrimeOrder
  V-Sorting samples
  V-Sorting samples time: 00:00:39
  Allocating rank array
  Ranking v-sort output
  Ranking v-sort output time: 00:00:21
  Invoking Larsson-Sadakane on ranks
  Invoking Larsson-Sadakane on ranks time: 00:00:37
  Sanity-checking and returning
Building samples
Reserving space for 172 sample suffixes
Generating random suffixes
QSorting 172 sample offsets, eliminating duplicates
QSorting sample offsets, eliminating duplicates time: 00:00:00
Multikey QSorting 172 samples
  (Using difference cover)
  Multikey QSorting samples time: 00:00:00
Calculating bucket sizes
Splitting and merging
  Splitting and merging time: 00:00:00
Split 22, merged 73; iterating...
Splitting and merging
  Splitting and merging time: 00:00:00
Split 8, merged 13; iterating...
Splitting and merging
  Splitting and merging time: 00:00:00
Split 3, merged 5; iterating...
Splitting and merging
  Splitting and merging time: 00:00:00
Split 2, merged 2; iterating...
Splitting and merging
  Splitting and merging time: 00:00:00
Split 1, merged 1; iterating...
Avg bucket size: 2.30677e+07 (target: 31087306)
Converting suffix-array elements to index image
Allocating ftab, absorbFtab
Entering GFM loop
Getting block 1 of 115
  Reserving size (31087307) for bucket 1
Getting block 2 of 115
  Reserving size (31087307) for bucket 2
  Calculating Z arrays for bucket 1
Getting block 3 of 115
  Reserving size (31087307) for bucket 3
  Calculating Z arrays for bucket 2
Getting block 4 of 115
  Reserving size (31087307) for bucket 4
Getting block 5 of 115
  Reserving size (31087307) for bucket 5
  Entering block accumulator loop for bucket 2:
  Entering block accumulator loop for bucket 1:
Getting block 6 of 115
  Reserving size (31087307) for bucket 6
  Calculating Z arrays for bucket 3
  Calculating Z arrays for bucket 4
Getting block 7 of 115
  Reserving size (31087307) for bucket 7
  Calculating Z arrays for bucket 5
Getting block 8 of 115
  Reserving size (31087307) for bucket 8
  Entering block accumulator loop for bucket 4:
Getting block 9 of 115
  Reserving size (31087307) for bucket 9
  Entering block accumulator loop for bucket 3:
  Calculating Z arrays for bucket 7
Getting block 10 of 115
  Calculating Z arrays for bucket 8
  Reserving size (31087307) for bucket 10
  Calculating Z arrays for bucket 9
  Entering block accumulator loop for bucket 5:
Getting block 11 of 115
  Reserving size (31087307) for bucket 11
  Entering block accumulator loop for bucket 8:
  Calculating Z arrays for bucket 10
  Entering block accumulator loop for bucket 9:
  Entering block accumulator loop for bucket 7:
  Calculating Z arrays for bucket 11
  Calculating Z arrays for bucket 6
  Entering block accumulator loop for bucket 10:
  Entering block accumulator loop for bucket 11:
  Entering block accumulator loop for bucket 6:
Getting block 12 of 115
  Reserving size (31087307) for bucket 12
  Calculating Z arrays for bucket 12
  Entering block accumulator loop for bucket 12:
Getting block 13 of 115
  Reserving size (31087307) for bucket 13
Getting block 14 of 115
  Reserving size (31087307) for bucket 14
  Calculating Z arrays for bucket 14
  Entering block accumulator loop for bucket 14:
Getting block 15 of 115
  Reserving size (31087307) for bucket 15
  Calculating Z arrays for bucket 15
  Entering block accumulator loop for bucket 15:
Getting block 16 of 115
  Reserving size (31087307) for bucket 16
  Calculating Z arrays for bucket 16
  Entering block accumulator loop for bucket 16:
  Calculating Z arrays for bucket 13
  Entering block accumulator loop for bucket 13:
  bucket 2: 10%
  bucket 3: 10%
  bucket 4: 10%
  bucket 6: 10%
  bucket 10: 10%
  bucket 5: 10%
  bucket 7: 10%
  bucket 12: 10%
  bucket 14: 10%
  bucket 8: 10%
  bucket 16: 10%
  bucket 4: 20%
  bucket 3: 20%
  bucket 6: 20%
  bucket 2: 20%
  .......
  .......
Wrote 1165549977 bytes to primary GFM file: genome.5.ht2
Wrote 675339278 bytes to secondary GFM file: genome.6.ht2
Re-opening _in5 and _in5 as input streams
Returning from HGFM constructor
Headers:
    len: 2652783500
    gbwtLen: 2652783501
    nodes: 2652783501
    sz: 663195875
    gbwtSz: 663195876
    lineRate: 6
    offRate: 4
    offMask: 0xfffffff0
    ftabChars: 10
    eftabLen: 0
    eftabSz: 0
    ftabLen: 1048577
    ftabSz: 4194308
    offsLen: 165798969
    offsSz: 663195876
    lineSz: 64
    sideSz: 64
    sideGbwtSz: 48
    sideGbwtLen: 192
    numSides: 13816581
    numLines: 13816581
    gbwtTotLen: 884261184
    gbwtTotSz: 884261184
    reverse: 0
    linearFM: Yes
Total time for call to driver() for forward index: 00:27:10
[jl19@spectre14 mm10]$ 
```
A number of files with .ht2 extension would be created. They are the index files
```
/scratch/bbash/jl19/Data_QC/mm10/genome.1.ht2
/scratch/bbash/jl19/Data_QC/mm10/genome.2.ht2
/scratch/bbash/jl19/Data_QC/mm10/genome.3.ht2
/scratch/bbash/jl19/Data_QC/mm10/genome.4.ht2
/scratch/bbash/jl19/Data_QC/mm10/genome.5.ht2
/scratch/bbash/jl19/Data_QC/mm10/genome.6.ht2
/scratch/bbash/jl19/Data_QC/mm10/genome.7.ht2
/scratch/bbash/jl19/Data_QC/mm10/genome.8.ht2
```
[The script](hisat2-creating-index-script.txt) for download MM10 and build and index the genome.fa.

[There are also prebuilt indexes on HISAT2's website ](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data)


For alignment, only mm10/genome is required.
```  
hisat2 -p 4 --dta -x mm10/genome -U trimmed/SRR7457551-trimmed.fastq.gz -S SRR7457551.sam;
hisat2 -p 4 --dta -x mm10/genome -U trimmed/SRR7457552-trimmed.fastq.gz -S SRR7457552.sam;
hisat2 -p 4 --dta -x mm10/genome -U trimmed/SRR7457555-trimmed.fastq.gz -S SRR7457555.sam;
hisat2 -p 4 --dta -x mm10/genome -U trimmed/SRR7457556-trimmed.fastq.gz -S SRR7457556.sam;
hisat2 -p 4 --dta -x mm10/genome -U trimmed/SRR7457559-trimmed.fastq.gz -S SRR7457559.sam;
hisat2 -p 4 --dta -x mm10/genome -U trimmed/SRR7457560-trimmed.fastq.gz -S SRR7457560.sam;
```  

output sam files (e.g. SRR7457551.sam) are located under /scratch/bbash/jl19/Data_QC/

Estimate run time is 20 mins for each run.


## Converting SAM to BAM using samtools "view"

``` 
module load samtools/1.15
``` 
``` 
samtools view --threads 8 -S -b SRR7457551.sam -o SRR7457551.bam;
samtools view --threads 8 -S -b SRR7457552.sam -o SRR7457552.bam;
samtools view --threads 8 -S -b SRR7457555.sam -o SRR7457555.bam;
samtools view --threads 8 -S -b SRR7457556.sam -o SRR7457556.bam;
samtools view --threads 8 -S -b SRR7457559.sam -o SRR7457559.bam;
samtools view --threads 8 -S -b SRR7457560.sam -o SRR7457560.bam;
``` 
## Sort the bam files
``` 
samtools sort SRR7457551.bam -o SRR7457551.sorted.bam;
samtools sort SRR7457552.bam -o SRR7457552.sorted.bam;
samtools sort SRR7457555.bam -o SRR7457555.sorted.bam;
samtools sort SRR7457556.bam -o SRR7457556.sorted.bam;
samtools sort SRR7457559.bam -o SRR7457559.sorted.bam;
samtools sort SRR7457560.bam -o SRR7457560.sorted.bam;
``` 
## Index sorted bam files
``` 
samtools index SRR7457551.sorted.bam;
samtools index SRR7457552.sorted.bam;
samtools index SRR7457555.sorted.bam;
samtools index SRR7457556.sorted.bam;
samtools index SRR7457559.sorted.bam;
samtools index SRR7457560.sorted.bam;
``` 
## Index files are generated

``` 
SRR7457551.sorted.bam.bai;
SRR7457552.sorted.bam.bai;
SRR7457555.sorted.bam.bai;
SRR7457556.sorted.bam.bai;
SRR7457559.sorted.bam.bai;
SRR7457560.sorted.bam.bai;
``` 
Bash script can be used to automate the procedure, the jobs can be submitted on Spectre.  This is beyond what this course can include.

## Check the BAM files

```
samtools view SRR7457551.bam  | head
``` 
``` 
[jl19@spectre14 Data_QC]$ samtools view SRR7457551.bam  | head
SRR7457551.2.1	16	chr15	100644892	60	74M	*	0	0GCTGGTTCTGGATTAGAAGAGGCGTCTGGGGAGGCCCTGGGGCTGCACAGCTGTGAAGAACCCTCTGGAGAGCC	??????????????????????????????????????????????????????????????????????????	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:74	YT:Z:UU	NH:i:1
SRR7457551.3.1	4	*	0	0	*	*	0	0	CCGTCCCGGATGGCCTTGGCGACAATGAACTCTGCATCTTCTGGGCTATCCAGCTGCAGCTTCTGGGAGATGTCGG	????????????????????????????????????????????????????????????????????????????	YT:Z:UU
SRR7457551.1.1	16	chr15	99409657	60	75M	*	0	0AAGGTGCCGGCCGCCCTGCCTGCCCGGTGCCTCCCACTGACCTGAGCTTACTAAGGGATGGCAGCAGCTCCCCAG	???????????????????????????????????????????????????????????????????????????	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:75	YT:Z:UU	NH:i:1
SRR7457551.5.1	4	*	0	0	*	*	0	0	GTCAGGGTTACCATGGCAACTGCATCATGGTTGGTAATTTTATATGTCCGACACCAAGAAATATGTCAGGGTTACC	????????????????????????????????????????????????????????????????????????????	YT:Z:UU
SRR7457551.6.1	0	chr9	50915472	60	73M	*	0	0GCTGTGCATCAGTCGAAGAACCTGCTCATTAAATTCTCCAATGGACGGCTCGTTTTCCTGTTCTTGTGGGTAG	?????????????????????????????????????????????????????????????????????????	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:73	YT:Z:UU	NH:i:1
SRR7457551.7.1	16	chr1	184731976	60	74M	*	0	0CCCCGACGCCGGCGTGCAGGATGTCCGCGATGCAGAAAGAGGGTTTCTTGACAGCGGCTGGGTCCAGGGCGAAG	??????????????????????????????????????????????????????????????????????????	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:74	YT:Z:UU	NH:i:1
SRR7457551.4.1	0	chr4	21776279	60	39M4023N33M	*	00	GAAGCTCCTGCGTGTGGAAGCTGCGGCCCGGCGGGCGGGTAAATAACAGATGCGGGTAAAAGATCCATCAAA????????????????????????????????????????????????????????????????????????	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:72	YT:Z:UU	XS:A:+	NH:i:1
SRR7457551.8.1	16	chr11	49683203	60	76M	*	0	0CCATGATGGACAGTAGCCGTATAACTGCCGGGTCGCATATTTATCACAAAGAACATTATAGCACATGACAGAAAAC	????????????????????????????????????????????????????????????????????????????	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:76	YT:Z:UU	NH:i:1
SRR7457551.9.1	0	chr5	120484110	60	76M	*	0	0GTGCTCCCCTCTGGGGTCTTGGGTCTCCCAGTGGATTGACTGGTGGAACCCACACATCCGTAAGTCAGGAGAGAGG	????????????????????????????????????????????????????????????????????????????	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:76	YT:Z:UU	NH:i:1
SRR7457551.12.1	16	chr19	4104891	60	74M2S	*	0	0	GACATCCGGGCACTGGAGGAGGAGACTGCACGCATGCTAGCGCAGCGTATGGCTAAGTGCAACACTGGCAGTGAAG	????????????????????????????????????????????????????????????????????????????	AS:i:-2	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:74	YT:Z:UU	NH:i:1
``` 

## View BAM files using IGV

The Integrative Genomics Viewer (IGC) is a high-performance, easy-to-use, intereactive tool for the visaul exploration of genomic data.

It can be downloaded as a desktop application or use as a web application.  It is also available to HPC

``` 
#load module
module load igv/2.3.82

Loading igv/2.3.82
  Loading requirement: java/1.8

#run igv
igv
``` 



[Go to Day1 Practical 3: Pseudoaligner and RNA-Seq Quantification Tool:kallisto](kallisto_alignment.md)

## [Back to Day 1 Practicals](rna-seq-wes-data-analysis-day1.md/#day-1-practicals)
