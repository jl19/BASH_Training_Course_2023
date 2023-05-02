# RNA-Seq Data Analysis

## [Day2](rna-seq-wes-data-analysis-day2.md)


> -  Quantification
> -  Differential Expression
> -  Visualize the output


## Quantification

Read counting is an important setp in quantifying gene expression levels from RNA-Seq data.  It involves counting the number of read that overlap with each gene/transcript after alignment step and quantifying their abundance.  RNA-Seq read counts often need to be normalized to account for differences in library size, sequencing depth as well as other technical factors.  Many software tools have been developed for quantification.

|Name|Description|
|----|----|
|[HTSeq](https://htseq.readthedocs.io/en/master/index.html) | HTSeq is a Python-based tool for counting reads mapped to features such as genes or exons. It can be used with a wide range of annotation formats and can handle single-end or paired-end reads.|
|[FeatureCounts](https://subread.sourceforge.net/featureCounts.html) | featureCounts is part of the Subread package. It can count reads mapped to genes, exons, or other features and can handle both single-end and paired-end reads.|
|[Salmon](https://combine-lab.github.io/salmon/)| Salmon is a tool to quantify the expression of transcripts of RNA-Seq data. It provides accurate expression estimates very quickly while using less memory.|
|[Kallisto](https://pachterlab.github.io/kallisto/about)|is a tool for quantifying abundances of transcripts from bulk and single-cell RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads. It is based on the novel idea of pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment. It is fast and memory-efficient and can handle large-scale RNA-Seq datasets.
|[StringTie](https://ccb.jhu.edu/software/stringtie/)|StringTie is a transcriptome assembly and quantification tool that can identify novel transcripts and splice variants. It can also quantify expression levels of known transcripts and is useful for gene fusions and other complex transcript structures.|
|[RSEM](https://deweylab.github.io/RSEM/README.html) |RSEM (RNA-Seq by Expectation-Maximization) is a tool for quantifying transcript expression that uses a statistical model based on the expectation-maximization algorithm. It can estimate expression levels of isoforms and genes and can handle both single-end and paired-end reads.|
|[Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)|Cufflinks is a popular tool for transcriptome assembly and quantification that can identify novel transcripts and splice variants. It can also quantify expression levels of known transcripts and is useful for gene fusions and other complex transcript structures.|


>  Here we use StringTie as an example for quantification.

>  [Practical 4a: Quantification using StringTie](practical-expression-quantification.md)



## Differential Expression Analysis 

Differentially expressed gene expression analysis (DE analysis) is a technique used to compare the gene/transcript expression levels of different samples or groups of samples. The goal of DE analysis is to identify genes that show significant changes in expression between the compared conditions.

DE analysis typically involves the following steps:

|Process|Description|
|----|----|
|Data preprocessing|Raw gene expression data is processed to remove systematic biases and normalize the data|
|Statistical analysis|Statistical tests are performed to identify genes that are differentially expressed between the conditions being compared. Commonly used tests include [t-tests](https://www.scribbr.com/statistics/t-test/), [ANOVA](https://www.scribbr.com/statistics/one-way-anova/), and the [Wilcoxon rank-sum test](https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_nonparametric/bs704_nonparametric4.html)|
|Multiple testing correction| Since thousands of genes are typically analyzed in a single experiment, statistical corrections are applied to adjust for the increased probability of false positives due to multiple comparisons|
|Functional analysis| Genes that are identified as differentially expressed are often analyzed to determine their biological functions and pathways|

DE analysis can provide insight into the molecular mechanisms underlying various biological processes and diseases. It is widely used in areas such as cancer research, drug discovery, and personalized medicine.

There are several software tools available for differential gene expression analysis. The choice of tool depends on various factors, such as the type and size of the data, the level of expertise of the user, and the availability of computational resources. Here are some commonly used tools for DE analysis:

|Tool Name|Description|
|----|----|
|[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) | A widely used [R/Bioconductor package](https://www.bioconductor.org/) for differential gene expression analysis. It is suitable for RNA-seq data and provides a robust method for normalization and variance stabilization|
|[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)| Another popular R/Bioconductor package for differential gene expression analysis. It is suitable for RNA-seq data and provides a flexible method for modeling gene expression data|
|[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)| A Bioconductor package that is widely used for microarray data analysis. It provides a comprehensive framework for preprocessing, normalization, and statistical analysis of microarray data|
|[limma voom](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html)| A Bioconductor package that provides a transformation method for RNA-seq data to allow for the application of linear models for statistical analysis|
|[NOISeq](https://bioconductor.org/packages/release/bioc/html/NOISeq.html)| A non-parametric method that is suitable for small sample sizes and can be used with both RNA-seq and microarray data|
|[DEGseq](https://bioconductor.org/packages/release/bioc/html/DEGseq.html)| A package that can be used with both RNA-seq and microarray data and provides a method for identifying differentially expressed genes based on negative binomial distribution|
|[Cuffdiff2](https://chipster.csc.fi/manual/cuffdiff2.html)| A package that is specifically designed for RNA-seq data and provides a method for identifying differentially expressed genes and isoforms|

> Differential Expression analysis using DESeq2

>  [Practical 5: Differential Expression analysis using DESeq2](analyzing-RNA-seq-data-with-DESeq2.md)

Compare sample groups: differential expression analysis


## Visualization of DE Outputs

Visualization of differential gene expression (DE) analysis output is an essential tool in the DE data analysis, it helps to understand the results and identify patterns in the data. Some commonly used visualization techniques for DE analysis output:

|Plot Name|Description|
|----|----|
|Volcano plots|A volcano plot is a scatter plot that shows the fold change of each gene on the x-axis and the statistical significance on the y-axis. Genes that are significantly differentially expressed are plotted above a threshold (e.g., p-value < 0.05) and are colored by their log-fold change|
|Heatmaps| Heatmaps are graphical representations of gene expression data that use color-coded rectangles to represent the expression level of each gene across different conditions or samples. Heatmaps can help to identify groups of genes that have similar expression patterns|
|MA plots| An MA plot is a scatter plot that shows the log-intensity ratio (M-value) on the y-axis and the average log-intensity (A-value) on the x-axis. It is a useful tool for visualizing the differences in gene expression between conditions|
|PCA plots| A PCA (principal component analysis) plot is a dimensionality reduction technique that can be used to visualize the overall patterns of gene expression across multiple samples or conditions. It can help to identify any clustering or separation of samples based on gene expression patterns|
|Pathway enrichment plots| Pathway enrichment analysis identifies the biological pathways that are significantly affected by differential gene expression. Pathway enrichment plots can help to visualize the pathways that are affected by differential gene expression and the genes that are involved in these pathways|

Visualization can be acheived using various software tools such as R/Bioconductor packages (e.g., ggplot2, pheatmap, PCAtools), Python packages (e.g., Matplotlib, Seaborn), and web-based tools (e.g., WEBGestalt GSEA, DAVID). 


>  [Practical 6: Visualization of DE Output](analyzing-RNA-seq-data-with-DESeq2-data-visualization.md)

## [Go to Day1 Data Analysis](rna-seq-wes-data-analysis-day1.md)
