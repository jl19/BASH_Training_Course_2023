
#STAR for Alignment

STAR  (Spliced Transcripts Alignment to a Reference). STAR is an aligner designed to specifically address many of the challenges of RNA-seq data mapping using a strategy to account for spliced alignments.


STAR has high accuracy but it is memory intensive.

```
module load star/2.7.1a 

```

Aligning reads using STAR is a two step process:

1. Create a genome index
2. Map reads to the genome

[STAR manual](/https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
##BWA (Burrows-Wheeler Aligner)

BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM.

Most commonly-used
Can do very long reads

Very quck alignment tool
Match  alignment positions with known gene positions


```
module load 
bwa/0.7.17 

```

BWA is top choice of aligments for downstream variant calling because it is accuracy.


Aligning reads with BWA-MEM

1. Creating BWA-MEM index

2. Aligment




##Somatic short varinat discovery (SNVs + Indels)

GATK Best Practices for data pre-processing

Call candidate variants: Mutect2

Annotate Variants:  Funcotator
Generate a list of mutation/variation

Highlight Nucleotide Variations, relative to the Referene Genome

Popular Visualization for Variation Calling Tools


IGV (Integrative Genomics Viewer)
GBrowse

JBrowse: Successor to GBrowse




Count how many rads each gene has




Q = − 10 log 10 P 20 1 in 100 99%
For example, if Phred assigns a Q score of 30 (Q30) to a base, this is
equivalent to the probability of an incorrect base call 1 in 1000 times
(Table 1). This means that the base call accuracy (i.e., the probability of
a correct base call) is 99.9%. A lower base call accuracy of 99% (Q20)
will have an incorrect base call probability of 1 in 100, meaning that
every 100 bp sequencing read will likely contain an error. When se-
quencing quality reaches Q30, virtually all of the reads will be perfect,
having zero errors and ambiguities. This is why Q30 is considered a
benchmark for quality in next-generation sequencing. 


Phred Quality Score          Probability of Incorrect Base Call    Base Call Accuracy
10                              1 in 10                             90%
20                              1 in 100                            99%
30                              1 in 1,000                          99.9%
40                              1 in 10,000                         99.99%
50                              1 in 100,000                        99.999%


### Pre-processing (FASTQC, Skewer, Trimmomatic)



Use Trimmomatic/Skewer to trim/remove poor quality bases/reads
Use Trimmomatic/Skewer to remove 3' adapter sequences from reads

| File Name   |File Extension | File Type  | Description                                      |                       
| ----------- | --------------|------------|--------------------------------------------------|
| `Fasta`     | .fasta, .fa   |sequences   |txt file for nucleotie or peptie sequencese       |
| `Fastq`     | .fastq, .fq   |read data   |txt file storing both sequence and its quality scores |
| `SAM`       | .sam          |short read alignment|sequence alignment file                   |
| `BAM`       |.bam           | binary SAM |
| `VCF`       | .vcf          | Variant information, Variant Call Format| txt file for variants calls incluing aNP, indels, CNV
|`GFF/GTF`    | .gff, .gtf    | annotation data | General Transfer/Feature Format annotation file


|Tool Name                                              | Description                                                   |
|--------------------|-----------------------|
|[Skewer](github.com/relipmoc/skewer)                   | A fast and sensitive trimmer for illumina paired-end sequences.|
|[Trimmomatic](https://github.com/usadellab/Trimmomatic)|Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.|
|[fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)|The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.
|
|[cutadapt](https://cutadapt.readthedocs.io/en/stable/)|Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.




[Practical 1: Pre-processing](practical-fastqc-2.md) 

Use [NoMachine NX Client to login to SPECTRE](NoMachine-login.md) 



###Align reads to a reference genome (HISAT2, BAM)
[Practical 2: Align reads to a reference genome](align-reads-to-reference-genome.md)


### 
[Practical 3: Quantification](practical-expression-quantification.md)


### 

[Practical 4: Differential Expression analysis (DESeq2, EdgeR, text files)](practical-differentail-expression-analysis.md)

Compare sample groups: differential expression analysis

Practicals:
