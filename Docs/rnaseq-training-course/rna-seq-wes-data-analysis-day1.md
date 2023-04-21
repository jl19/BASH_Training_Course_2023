# RNA-Seq Data Analysis


## [Day1](rna-seq-wes-data-analysis-day1.md)
:memo:
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

### Advantages of NGS over traditional methods:

* Higher sensitivity to detect low-frequency variants
* Higher genomic coverage
* To be able to sequence hundreds to thousands of genes or gene regions simultaneously 
* Higher sample throughput

## RNA-Seq Data Analysis: Steps, Tools and File formats


RNA-Seq (RNA-sequencing) refers to analyze the transcriptome of biological sample using next-generation sequencing (NGS) technologies. It can  be used to quantify the abundance of RNA as well as identify alternative splicing events, gene fusion etc. It can also be applied in detect non-coding RNAs e.g. microRNAs and long non-coding RNA.
### RNA-Seq Workflow

![RNA-Seq Work Flow ](https://jl19.github.io/BASH_Training_Course_2023/Docs/assets/RNA-Seq_work_flow.jpeg)


### Perform QC for better downstream analysis (FastQC, FASTQ files)

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

* To identify potential issues before data analysis

*  Base calling accuracy, measured by the Phred quality score (Q score), is the most common metric used to assess the accuracy of a sequencing platform.  It indicates the probablity that a given base is called incorrectly by the sequencer.

>:bulb:  Q = − 10 log 10(e)

Where e is the estimated probability of the base call being wrong.

* Higher Q scores indicate a smaller prbability of error.
* Lower Q scores can results in significant portion of the reads being unusable.  They may also lead to increase false-positive variant calls, resulting in inaccurate conclusions.

For example, if Phred assigns a Q score of 30 (Q30) to a base, this is equivalent to the probability of an incorrect base call 1 in 1000 times
(see Table below). 

This means that the base call accuracy (i.e., the probability of a correct base call) is 99.9%. 

|Phred Quality Score          |Probability of Incorrect Base Call |   Base Call Accuracy|
| ----------- |---------------|--------------------------------------------------|
|10                            |  1 in 10                        |       90%|
|20                            |  1 in 100                       |       99%|
|30                            |  1 in 1,000                     |     99.9%|
|40                            |  1 in 10,000                    |     99.99%|
|50                            |  1 in 100,000                   |    99.999%|

Q20 (99%) will have an incorrect base call probability of 1 in 100, it represents every 100 bp sequencing read will likely contain an error. When sequencing quality reaches Q30, virtually all of the reads will be perfect, having zero errors and ambiguities. Q30 is considered a benchmark for quality in next-generation sequencing. 

### Pre-processing (FASTQC, Skewer, Trimmomatic)

#### Why trim RNA-Seq data 

Read trimming tools have been developed to remove adapter sequences corresponding to the library adapters present in the FASTQ files and bases with low sequencing quality from sequencing reads such as RNA-seq reads, in order to help read aligners to achieve a better read mapping result.

Use Trimmomatic/Skewer to trim/remove poor quality bases/reads

Use Trimmomatic/Skewer to remove 3' adapter sequences from reads

>:bulb: Pre-processing Tools:

|Tool Name  | Description   |
|--------------------|-----------------------|
|[Skewer](https://github.com/relipmoc/skewer)                   | A fast and sensitive trimmer for illumina paired-end sequences.|
|[Trimmomatic](https://github.com/usadellab/Trimmomatic)|Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.|
|[fastx_toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)|The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.
|[cutadapt](https://cutadapt.readthedocs.io/en/stable/)|Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

Here are the common file formats in NGS:


| File Name   |File Extension | File Type  | Description                                      |                       
| ----------- | --------------|------------|--------------------------------------------------|
| [Fasta](FASTA-FORMAT.md)| .fasta, .fa   |sequences   |txt file for nucleotie or peptie sequences       |
| [Fastq](https://learn.gencore.bio.nyu.edu/ngs-file-formats/fastq-format/)     | .fastq, .fq   |read data   |txt file storing both sequence and its quality scores |
| [SAM](https://learn.gencore.bio.nyu.edu/ngs-file-formats/sambam-format/)       | .sam          |short read alignment|sequence alignment file         |
| [BAM](https://learn.gencore.bio.nyu.edu/ngs-file-formats/sambam-format/)      |.bam           | binary SAM |binary which means that it is significantly smaller than the SAM files and significantly faster to read| 
|[VCF](https://learn.gencore.bio.nyu.edu/ngs-file-formats/vcf-format/)| .vcf| Variant information, Variant Call Format| txt file for variants calls incluing SNP, indels, CNV
|[GFF/GTF](https://learn.gencore.bio.nyu.edu/ngs-file-formats/gff3-format/)| .gff, .gtf| annotation data | General Transfer/Feature Format annotation file



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

>:bulb: Commonly Used Alignment Tools

|Alignment Tool Name| Description|
| ----------- | ------------------------------------ |
|[Burrows-Wheeler Aligner (BWA)](https://github.com/lh3/bwa) |a widely used algorithm that utilizes the Burrows-Wheeler transform to index the reference genome and align the reads|
|[TopHat](https://ccb.jhu.edu/software/tophat/index.shtml) |a splice-aware alignment algorithm that can detect spliced alignments and map reads across splice junctions|
|[STAR](https://github.com/alexdobin/STAR) |a fast and accurate alignment algorithm that can align reads to both the genome and transcriptome simultaneously|
|[HISAT2](http://daehwankimlab.github.io/hisat2/) |a newer splice-aware alignment algorithm that uses a hierarchical indexing strategy to improve speed and accuracy|
|[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) | a fast and memory-efficient algorithm that aligns reads to the reference genome using a seed-and-extend approach|
|[Subread](https://subread.sourceforge.net/) | an alignment algorithm that uses an efficient seed-and-vote strategy to align reads to the reference genome|

Each of these algorithms has its own strengths and weaknesses, and the choice of algorithm depends on factors such as the size of the reference genome, the sequencing technology used, and the research question being addressed.

#### Pseudoaligment


Pseudoaligment is to align short reads generated by NGS to a reference genome without explicitly perfoming full read alignment. It uses an indexing strategy that can quickly dtermine which genomic regions or transcripts each read aligns to in order to achive faster analysis.

|Pseudoalignment Tool Name| Description|
| ----------- | ------------------------------------ |
|[Kallisto](https://pachterlab.github.io/kallisto/about) | Uses a k-mer based approach to build an index of transcripts and assing reads to the transcripts base on shared k-mers|
|[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) | It uses a lightweight alignment strategy to map the reads to the transcripts in the reference genome or transcriptome, it quantifies gene expression levels using a statistical model that takes into account the fragment-level bias in the sequencing data. It is useful when comparing multiple samples simultaneously|

Pseudoaligner is suitable for larger dataset due to its speed but it is not suitable for detecting of novel transcripts and variant calling.

## RNA-Seq Applications

|Application| Description|
| ----------- | ------------------------------------ |
|Gene Differential Expression|used to compare gene expression levels between different conditions, such as healthy versus diseased tissue, or treated versus untreated cells in order to explore the molecular mechanisms underlying different biological processes|
|Transcriptome annotation| To identify and annotate novel transcripts, including alternative splicing events and non-coding RNAs in order to understand the complexity of the transcriptome and its functional implications|
|Functional genomics|To identify genes involved in specific biological processes, such as stress response or cell differentiation. This can help researchers identify potential drug targets or therapeutic strategies for specific diseases|
|Biomarker discovery| To identify biomarkers, or molecular signatures, associated with specific diseases or conditions. This is useful for  developing diagnostic tools and identifying potential targets for personalized therapies|



---------------------------

## Day 1 Practicals

---------------------------------
Please following the instruction from Research Computing of University of Leicester to install the NoMachine NX Client [NoMachine NX Client to login to SPECTRE](NoMachine-login.md) 

### Pre-processing

> :boom: [Practical 1: Pre-processing](pre-processing-rna.md) 

In this practical you will learn :

* To import/copy, view and check the quality of the sequencing data
* To trim the raw reads using a trimming tool

-------------------------
### Alignment (HISAT2, STAR, BWA, Tophat, Kallisto, Salmon)

> Alignment to a reference genome using HISAT2

> :boom: [Practical 2: Align reads to a reference genome using HISAT2](align-reads-to-reference-genome.md)

In this practical you will learn to align reads to a reference genome using HISAT2, an alignment tool.

--------------

> Alignment to transcriptome using pseudoaligner: Kallisto 

> :boom: [Practical 3: Pseudoaligner and RNA-Seq Quantification Tool:kallisto](kallisto_alignment.md)

In this practical you will learn to use a pseudoalinger tool Kallistol to align and quantify.

----------------------

## [Go to Day2 Data Analysis](rna-seq-wes-data-analysis-day2.md)
