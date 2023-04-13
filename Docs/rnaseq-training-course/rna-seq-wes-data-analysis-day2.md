# RNA-Seq Data Analysis

## [Day2](rna-seq-wes-data-analysis-day2.md)


> -  Quantification
> -  Differential Expression
> -  Visualize the output


## Quantification

Read counting is an important setp in quantifying gene expression levels from RNA-Seq data.  It involves counting the number of read that overlap with each gene/transcript after alignment step and quantifying their abundance.  RNA-Seq read counts often need to be normalized to account for differences in library size, sequencing depth as well as other technical factors.  Many software tools have been developed for quantification.


HTSeq
FeatureCounts
Salmon
Kallisto
StringTie
RSEM
Cufflinks

Genes/Transcripts counting tools can be used to obtain abundance of genes/transcripts 

* StringTie 

* Kallisto

> :boom: [Practical 4: Quantification using StringTie](practical-expression-quantification.md)

## Differential Expression Analysis using DESeq2

Normalization and statistical testing to identify differentially expressed genes.

> :boom: [Practical 5: Differential Expression analysis using DESeq2](analyzing-RNA-seq-data-with-DESeq2.md)

Compare sample groups: differential expression analysis


## Visualization of the Outputs

> :boom: [Practical 6: Visualization](analyzing-RNA-seq-data-with-DESeq2-data-visualization.md)

## [Go to Day1 Data Analysis](rna-seq-wes-data-analysis-day1.md)
