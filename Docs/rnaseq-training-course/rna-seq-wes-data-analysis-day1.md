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

## RNA-seq Data Analysis: Steps, Tools and File formats


RNA-seq (RNA-sequencing) is a technique that can examine the quantity and sequences of RNA in a sample using next-generation sequencing (NGS). It analyzes the transcriptome, indicating which of the genes encoded in our DNA are turned on or off and to what extent. 

### RNA-Seq Workflow

![RNA-Seq Work Flow ](https://jl19.github.io/BASH_Training_Course_2023/Docs/assets/RNA-Seq_work_flow.jpeg)


### Perform QC for better downstream analysis (FastQC, FASTQ files)

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

* To identify potential issues before data analysis

*  Base calling accuracy, measured by the Phred quality score (Q score), is the most common metric used to assess the accracy of a sequencing platform.  It indicates the probablity that a given base is called incorrectly by the sequencer.



>  Q = − 10 log 10(e)

Where e is the estimated probability of the base call being wrong.

* Higher Q scores indicate a samller prbability of error.
* Lower Q scores can results in significant portion of the reads being unusable.  They may also lead to increase false-positive variant calls, resulting in inaccurate conclusions.

For example, if Phred assigns a Q score of 30 (Q30) to a base, this is
equivalent to the probability of an incorrect base call 1 in 1000 times
(see Table below). 

This means that the base call accuracy (i.e., the probability of
a correct base call) is 99.9%. 

|Phred Quality Score          |Probability of Incorrect Base Call |   Base Call Accuracy|
| ----------- |---------------|--------------------------------------------------|
|10                            |  1 in 10                        |       90%|
|20                            |  1 in 100                       |       99%|
|30                            |  1 in 1,000                     |     99.9%|
|40                            |  1 in 10,000                    |     99.99%|
|50                            |  1 in 100,000                   |    99.999%|

Q20 (99%) will have an incorrect base call probability of 1 in 100, it represents
every 100 bp sequencing read will likely contain an error. When sequencing quality reaches Q30, virtually all of the reads will be perfect, having zero errors and ambiguities. Q30 is considered a
benchmark for quality in next-generation sequencing. 

### Pre-processing (FASTQC, Skewer, Trimmomatic)

#### Why trim RNA-Seq data 

Read trimming tools have been developed to remove adapter sequences corresponding to the library adapters present in the FASTQ files and bases with low sequencing quality from sequencing reads such as RNA-seq reads, in order to help read aligners to achieve a better read mapping result.

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


### Alignment

RNA-Seq alignment is the process of aligning the short reads genearted from RNA sequencing to a reference genome or transcriptome.  It is used to identify the location of each read and determine the expression level of genes/transcripts.  There are various alignment tools and algorithms used for the alignment.  

* Indetify the best match/matches
* Read length
* Sequencing errors
* Genomic variation

The output of RNA-Seq alignment is a file that lists the location and quality of each read which can be used for downstream analysis e.g. gene expression quatification and differential gene expression analysis etc.

There are several types of alignment algorithms used for RNASeq alignment, including:

Burrows-Wheeler Aligner (BWA): a widely used algorithm that utilizes the Burrows-Wheeler transform to index the reference genome and align the reads.

TopHat: a splice-aware alignment algorithm that can detect spliced alignments and map reads across splice junctions.

STAR: a fast and accurate alignment algorithm that can align reads to both the genome and transcriptome simultaneously.

HISAT2: a newer splice-aware alignment algorithm that uses a hierarchical indexing strategy to improve speed and accuracy.

Bowtie2: a fast and memory-efficient algorithm that aligns reads to the reference genome using a seed-and-extend approach.

Subread: an alignment algorithm that uses an efficient seed-and-vote strategy to align reads to the reference genome.

Each of these algorithms has its own strengths and weaknesses, and the choice of algorithm depends on factors such as the size of the reference genome, the sequencing technology used, and the research question being addressed.

#### Pseudoaligment

Pseudoaligment is align short reads  generated by NGS to a reference genome without explicitly perfoming full read alignment. It uses an indexing strategy that can quickly dtermine which genomic regions or transcripts each read aligns to in order to achive faster analysis.

Kallisto: Uses a k-mer based approach to build an index of transcripts and assing reads to the transcripts base on shared k-mers.

Salmon: It uses a lightweight alignment strategy to map the reads to the transcripts in the reference genome or transcriptome, it quantifies gene expression levels using a statistical model that takes into account the fragment-level bias in the sequencing data. It is useful when comparing multiple samples simultaneously.

RapMap: It also uses a k-mer based indexing strategy to quickly identify potential transcript matches for each read. In addition, it incorporateds a vovel prbablilistic model that takes into account the uncertainty in the mapping process and the likelihood of observing each read given the transcripts it could  potentially originate from.

Pseudoaligner is suitable for larger dataset due to its speed but it is not suitable for detecting of novel transcripts and variant calling.

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
Please following the instruction from Research Computing of University of Leicester to install the NoMachine NX Client [NoMachine NX Client to login to SPECTRE](NoMachine-login.md) 

### Pre-processing

[Practical 1: Pre-processing](pre-processing-rna.md) 

In this practical you will learn :

* To import/copy, view and check the quality of the sequencing data
* To trim the raw reads using a trimming tool

-------------------------
### Alignment (HISAT2, STAR, BWA, Tophat, Kallisto, Salmon)

> Alignment to a reference genome using HISAT2

> [Practical 2: Align reads to a reference genome using HISAT2](align-reads-to-reference-genome.md)

In this practical you will learn to align reads to a reference genome using HISAT2, an alignment tool.

--------------

> Alignment to transcriptome using pseudoaligner: Kallisto 

> [Practical 3: Pseudoaligner and RNA-Seq Quantification Tool:kallisto](kallisto_alignment.md)

In this practical you will learn to use a pseudoalinger tool Kallistol to align and quantify.

----------------------


## [Go to Day2 Data Analysis](rna-seq-wes-data-analysis-day2.md)
