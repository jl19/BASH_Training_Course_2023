# RNA-Seq Data Analysis

## [Day2](rna-seq-wes-data-analysis-day2.md)


> -  Quantification
> -  Differential Expression
> -  Visualize the output


## Quantification

Read counting is an important setp in quantifying gene expression levels from RNA-Seq data.  It involves counting the number of read that overlap with each gene/transcript after alignment step and quantifying their abundance.  RNA-Seq read counts often need to be normalized to account for differences in library size, sequencing depth as well as other technical factors.  Many software tools have been developed for quantification.

|Name|Description|
|----|----|
|[HTSeq](https://htseq.readthedocs.io/en/master/index.html)HTSeq is a Python-based tool for counting reads mapped to features such as genes or exons. It can be used with a wide range of annotation formats and can handle single-end or paired-end reads.|
|FeatureCounts|featureCounts is part of the Subread package. It can count reads mapped to genes, exons, or other features and can handle both single-end and paired-end reads.|
|Salmon|Salmon is a lightweight alignment-free tool that quantifies transcript expression using k-mer-based mapping. It can accurately quantify transcript expression even for lowly expressed transcripts and can handle reads from multiple sequencing platforms.|
|Kallisto|Kallisto is another alignment-free tool that quantifies transcript expression using pseudoalignment. It is fast and memory-efficient and can handle large-scale RNA-Seq datasets.|
|StringTie|StringTie is a transcriptome assembly and quantification tool that can identify novel transcripts and splice variants. It can also quantify expression levels of known transcripts and is useful for gene fusions and other complex transcript structures.|
|RSEM|RSEM (RNA-Seq by Expectation-Maximization) is a tool for quantifying transcript expression that uses a statistical model based on the expectation-maximization algorithm. It can estimate expression levels of isoforms and genes and can handle both single-end and paired-end reads.|
|Cufflinks|Cufflinks is a popular tool for transcriptome assembly and quantification that can identify novel transcripts and splice variants. It can also quantify expression levels of known transcripts and is useful for gene fusions and other complex transcript structures.|


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
