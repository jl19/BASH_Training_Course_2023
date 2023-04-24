# DNA Library Preparation
A DNA library is a collection of DNA fragments that are specially prepared for sequencing. We are not going to get into the gory details of how libraries are prepared. All you need to know is that the end result of library prep is a bunch of DNA fragments that are made up of DNA from the organism you’re working with flanked by artificial DNA, which is added by you during library prep, which allows your DNA to stick to a solid surface and be sequenced on a sequencing machine. These artificial pieces of DNA are called adapters. 

## Single-End Read


![Single-End Read](https://jl19.github.io/BASH_Training_Course_2023/Docs/assets/Single-End-Read.jpg)

Single-end reads refer to a type of sequencing where only one end of a DNA fragment is sequenced. It reads in one direction from 5' prime end to 3' prime end.

Single-end sequencing has advantages in terms of cost, simplicity, and data analysis. It requires less sequencing depth and is easier to analyze, but may not provide as much information about the sequence as paired-end sequencing. 


The most basic adapter layout allows for single end sequencing. Elements 1 and 2 are the terminal parts of the adapter and are what bind to oligos, small pieces of DNA, present on a solid surface. These elements are also necessary for fragment amplification. 

The insert, the bit of DNA from our organism of interest, from element 1 we have a priming segment for read 1. The priming segment is recognized by a sequencing primer, which initiates the sequencing by synthesis process. 

Elements 1 and 2, as well as the Read 1 sequencing primer, are required by most sequencing technologies. 



## Paired-End Read

Paired-end (PE) sequencing generates reads from both ends of the fragment. The sequence first carried out from 5' prime to 3' prime direction, after turn around chemistry, it the sequence from 3' prime to 5' prime direction.  Two reads are read towards the middle of the insert.  Read 1 and read 2 normally are the same length.  Paired-end sequencing provides more information about the sequence and allows for the detection of structural variations in the genome.

![Paired-End Read](https://jl19.github.io/BASH_Training_Course_2023/Docs/assets/Paired-End-Read.jpg)


For paired-end read, the adapter is the same as the single end adapter, except it includes a second priming sequence for sequencing the other end of the read. This second priming sequence can also be used to read an index. An index allows for the pooling (mixing) of libraries from multiple samples together. Pooling enables the sequencing of multiple samples on the same lane of a sequencer, a practice often referred to as multiplexing.  The indices are used for demultiplexing, which is just the sorting out of reads based on their index sequence after sequencing. The index is read in a second sequencing step using an index sequencing primer. Alternatively, the index sequence can be included in-line as part of the first read. When an index is included in the first read, it is often referred to as a barcode. This is depicted in the third adapter layout diagram. Barcodes are read as the beginning part of the read and need to cut out after demultiplexing.

Paired-end sequencing provides several advantages over single-end sequencing. By sequencing both ends of the fragment, it allows for more accurate read alignment and improves the detection of structural variations such as insertions, deletions, and inversions. Paired-end sequencing can also increase the accuracy of base calling, which is the process of determining the nucleotide sequence of each read.
In addition to these advantages, paired-end sequencing also enables the identification of gene fusions and alternative splicing events, which can be useful for studying gene expression and identifying disease-causing mutations.

However, paired-end sequencing is typically more expensive and generates more data than single-end sequencing, which may require more extensive computational resources for data analysis.

## Single Reads and Paired-End Reads Application

|Application| Type| Description|
|----|----|----|
|SNP Detection| Single Read or Paired-End Read|Coverage depth is key|
|Indel or Stucture Variant Detection| Paired-End Read|Analysis methods are base on PE data|
|De Novo Genome or Transcriptome Assembly|Paired-End Read| PE info is used in assembly process|
|RNA-Seq (Expression)| Single or Paired-End| PE needed for identification of novel transcripts and gene structure characterization|
|Small RNA Differential Expression| Single Read|PE will result in high overlap|
