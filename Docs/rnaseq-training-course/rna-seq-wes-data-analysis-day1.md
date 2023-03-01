# RNA-seq (WES) Data Analysis


## [Day1](rna-seq-wes-data-analysis-day1.md)

> -  Understand the basics of the NGS technologies     
> -  RNA-seq Data Analysis: Steps, Tools and File formats
> -  RNA-seq applications
> -  Perform quality control for better downstream analysis
> -  Pre-processing
> -  Align reads to a reference genome


## Understand the basics of NGS technologies

Next-generation sequencing (NGS) is used to determine the order of nucleotides in entire genomes or targeted regions of DNA or RNA.

### The typical NGS process involves:

* Fragmenting DNA/RNA into multiple pieces
* Adding adapters
* Sequencing the libraries
* Reassembling them to form genomic sequence

### Advantages of NGS:

* Higher sensitivity to detect low-frequency variants
* Higher genomic coverage
* To be able to sequence hundreds to thousands of genes or gene regions simultaneously 
* Higher sample throughput

##RNA-seq Data Analysis: Steps, Tools and File formats


RNA-seq (RNA-sequencing) is a technique that can examine the quantity and sequences of RNA in a sample using next-generation sequencing (NGS). It analyzes the transcriptome, indicating which of the genes encoded in our DNA are turned on or off and to what extent. 

### RNA-Seq Workflow

![RNA-Seq Work Flow ](/Docs/assets/RNA-Seq_work_flow.jpeg)


### Perform QC for better downstream analysis (FastQC, FASTQ files)

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

* To identify potential issues before data analysis

*  Base calling accuracy, measured by the Phred quality score (Q score), is the most common metric used to assess the accracy of a sequencing platform.  It indicates the probablity that a given base is called incorrectly by the sequencer.



>  Q = − 10 log 10 P 20 1 in 100 99%

For example, if Phred assigns a Q score of 30 (Q30) to a base, this is
equivalent to the probability of an incorrect base call 1 in 1000 times
(see Table below). 

This means that the base call accuracy (i.e., the probability of
a correct base call) is 99.9%. 

|Phred Quality Score          |Probability of Incorrect Base Call |   Base Call Accuracy|
| ----------- |---------------|--------------------------------------------------|
|10                            |  1 in 10                       |      90%|
|20                            |  1 in 100                       |     99%|
|30                            |  1 in 1,000                     |     99.9%|
|40                            |  1 in 10,000                     |    99.99%|
|50                            |  1 in 100,000                    |    99.999%|

Q20 (99%) will have an incorrect base call probability of 1 in 100, it represents
every 100 bp sequencing read will likely contain an error. When sequencing quality reaches Q30, virtually all of the reads will be perfect,
having zero errors and ambiguities. Q30 is considered a
benchmark for quality in next-generation sequencing. 

### Pre-processing (FASTQC, Skewer, Trimmomatic)

Use Trimmomatic/Skewer to trim/remove poor quality bases/reads

Use Trimmomatic/Skewer to remove 3' adapter sequences from reads

> Pre-processing Tools:

|Tool Name                                              | Description                                                   |
|--------------------|-----------------------|
|[Skewer](https://github.com/relipmoc/skewer)                   | A fast and sensitive trimmer for illumina paired-end sequences.|
|[Trimmomatic](https://github.com/usadellab/Trimmomatic)|Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.|
|[fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)|The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.
|[cutadapt](https://cutadapt.readthedocs.io/en/stable/)|Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

File Formats:


| File Name   |File Extension | File Type  | Description                                      |                       
| ----------- | --------------|------------|--------------------------------------------------|
| [Fasta](https://www.bioinformatics.nl/tools/crab_fasta.html)    | .fasta, .fa   |sequences   |txt file for nucleotie or peptie sequencese       |
| [Fastq](https://en.wikipedia.org/wiki/FASTQ_format)     | .fastq, .fq   |read data   |txt file storing both sequence and its quality scores |
| [SAM](https://en.wikipedia.org/wiki/SAM_(file_format))       | .sam          |short read alignment|sequence alignment file                   |
| [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map)      |.bam           | binary SAM |
| [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf)       | .vcf          | Variant information, Variant Call Format| txt file for variants calls incluing aNP, indels, CNV
|[GFF/GTF](https://www.ensembl.org/info/website/upload/gff.html)    | .gff, .gtf    | annotation data | General Transfer/Feature Format annotation file


| File      | Description                          |
| ----------- | ------------------------------------ |
| `Single-end`       | Each  read is a single sequence from one end of a DNA fragment (single fastq file). The fragment is usually 200-800bp long, with the amount being read can be chosen between 50 and 300 bp  |
| `Paired-end`       | Each read is two sequences (a pair) from each end of the same DNA fragment. The distance between the reads on the original genome sequence is equal to the length of the DNA fragment that was sequenced, usually 200-800 bp |
| `Mate-pair`    | Each read is two sequences from each end of the same DNA fragment, but the distance between the reads on the original genome sequence is much longer, e.g. 3000-10000 bp


## RNA-Seq Applications

* Differential Expression
* Isoform switching
* New genes and transcripts
* New Transcriptiomes
* Variants Detection
* Allele-specific expression

---------------------------

## Day 1 Practicals

---------------------------------

### Pre-processing

[Practical 1: Pre-processing](pre-processing-rna.md) 

Use [NoMachine NX Client to login to SPECTRE](NoMachine-login.md) 


-------------------------
### Alignment (HISAT2, STAR, BWA, Kallisto)

> Alignment to a reference genome using HISAT2

> [Practical 2: Align reads to a reference genome using HISAT2](align-reads-to-reference-genome.md)

--------------

> Alignment to transcriptome using pseudoaligner: Kallisto 

> [Practical 3: Pseudoaligner and RNA-Seq Quantification Tool:kallisto](kallisto_alignment.md)

----------------------


## [Go to Day2 Data Analysis](/rnaseq-training-course/rna-seq-wes-data-analysis-day2/)
