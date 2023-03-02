# Analyzing RNA-seq data with DESeq2 using R

### Install required R packages
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
install.packages("reshape2")
library("DESeq2")
library("ggplot2")
library("reshape2")
```

### Loading data from previous step tximport output

```
##Construct a DESeqDataSEt from the txi object created in txtimport.R and sample information in samples.
#txi is matrices which contain the transcript-level abundance, estimated counts and transcript lengthes etc.


dds_txi<- DESeqDataSetFromTximport(txi,
                              colData = samples,
                              design = ~ condition)
#using counts and average transcript lengths from tximport
```
### Differential expression analysis

```
dds_txi <- DESeq(dds_txi)

# estimating size factors
# using 'avgTxLength' from assays(dds), correcting for library size
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
```
### Result for Navie vs 24h

```
# Results
res <- results(dds_txi)

res

# log2 fold change (MLE): condition Naive vs 24h 
# Wald test p-value: condition Naive vs 24h 
# DataFrame with 48321 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   0610006L08Rik   0.00000             NA        NA        NA          NA          NA
# 0610007P14Rik 447.24265      -1.555441  0.226567 -6.865263 6.63688e-12 1.99161e-10
# 0610009B22Rik 211.76994      -0.203472  0.242080 -0.840518 4.00618e-01 5.94800e-01
# 0610009E02Rik   1.15551       0.471663  4.453707  0.105904 9.15659e-01          NA
# 0610009L18Rik  19.00117       0.634050  0.746600  0.849250 3.95742e-01 5.90287e-01
# ...                 ...            ...       ...       ...         ...         ...
# Zyg11a             0.00             NA        NA        NA          NA          NA
# Zyg11b          1242.07       0.344672  0.135832   2.53748 1.11654e-02 4.09594e-02
# Zyx             5438.37       0.897308  0.127010   7.06486 1.60778e-12 5.24942e-11
# Zzef1           2635.28       0.178422  0.113657   1.56983 1.16455e-01 2.54259e-01
# Zzz3            1326.27       0.559339  0.154486   3.62064 2.93873e-04 1.92334e-03
```

### Compare two conditions using names
```
#Note that we could have specified the coefficient or contrast we want to build a results table for, 
#using either of the following equivalent commands:

#check existing names of the object dds_txi

resultsNames(dds_txi)
#[1] "Intercept"              "condition_2h_vs_24h"    "condition_Naive_vs_24h"

#Compare two conditions using name above

res_Naive_vs_24h <- results(dds_txi, name="condition_Naive_vs_24h")

#Compare two conditions using contrast
#Contrast is the factor set in tximport.R which are Navie, 2h and 24h

res_Naive_vs_2h <- results(dds_txi, contrast=c("condition","Naive","2h"))
```
### Result for Navie vs 2h

```
# log2 fold change (MLE): condition Naive vs 2h 
# Wald test p-value: condition Naive vs 2h 
# DataFrame with 48321 rows and 6 columns
# baseMean log2FoldChange     lfcSE       stat    pvalue      padj
# <numeric>      <numeric> <numeric>  <numeric> <numeric> <numeric>
#   0610006L08Rik   0.00000             NA        NA         NA        NA        NA
# 0610007P14Rik 447.24265       0.289074  0.233415   1.238454 0.2155478  0.503634
# 0610009B22Rik 211.76994       0.568372  0.248935   2.283212 0.0224179  0.107668
# 0610009E02Rik   1.15551      -3.993491  4.248810  -0.939908 0.3472647        NA
# 0610009L18Rik  19.00117       0.416083  0.746494   0.557383 0.5772655  0.838870
# ...                 ...            ...       ...        ...       ...       ...
# Zyg11a             0.00             NA        NA         NA        NA        NA
# Zyg11b          1242.07     -0.1872923  0.135213 -1.3851664  0.166002  0.431721
# Zyx             5438.37     -0.0872193  0.126580 -0.6890428  0.490796  0.781206
# Zzef1           2635.28     -0.0630209  0.113731 -0.5541236  0.579494  0.840554
# Zzz3            1326.27      0.0119787  0.154244  0.0776606  0.938098  0.986018

```

### Result for 24 vs 2h

```
res_24h_vs_2h <- results(dds_txi, contrast = c("condition", "24h", "2h"))

```

### Result for 24 vs Naive

```
res_24h_vs_Naive <- results(dds_txi, contrast=c("condition","24h", "Naive"))

```

### Result for 2h vs Navie

```
res_2h_vs_Naive <- results(dds_txi, contrast=c("condition","2h", "Naive"))
```


### Log fold change shrinkage for visualization and ranking

```
#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. 
#To shrink the LFC, we pass the dds object to the function lfcShrink. Below we specify to use the 
#apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the
#previous estimator.

#We provide the dds object and the name or number of the coefficient we want to shrink,
#where the number refers to the order of the coefficient as it appears in resultsNames(dds_txi)
```

### List Results Names

```
#using resultsNames to list the names of the estimated effects (coefficents) of the model
resultsNames(dds_txi)

## [1] "Intercept"                     "condition_2h_vs_24h"    "condition_Naive_vs_24h"

resNaive_vs_24h_LFC <- lfcShrink(dds_txi, coef="condition_Naive_vs_24h", type="apeglm")
resNaive_vs_24h_LFC

# log2 fold change (MAP): condition Naive vs 24h 
# Wald test p-value: condition Naive vs 24h 
# DataFrame with 48321 rows and 5 columns
# baseMean log2FoldChange     lfcSE      pvalue        padj
# <numeric>      <numeric> <numeric>   <numeric>   <numeric>
#   0610006L08Rik   0.00000             NA        NA          NA          NA
# 0610007P14Rik 447.24265     -1.4960322  0.230078 6.63688e-12 1.99161e-10
# 0610009B22Rik 211.76994     -0.1492869  0.211503 4.00618e-01 5.94800e-01
# 0610009E02Rik   1.15551     -0.0019808  0.389710 9.15659e-01          NA
# 0610009L18Rik  19.00117      0.1430094  0.373807 3.95742e-01 5.90287e-01
# ...                 ...            ...       ...         ...         ...
# Zyg11a             0.00             NA        NA          NA          NA
# Zyg11b          1242.07       0.314919  0.132775 1.11654e-02 4.09594e-02
# Zyx             5438.37       0.868430  0.127829 1.60778e-12 5.24942e-11
# Zzef1           2635.28       0.165588  0.110122 1.16455e-01 2.54259e-01
# Zzz3            1326.27       0.516176  0.154057 2.93873e-04 1.92334e-03
```

### Ordering results using p-values and adjusted p-values  

```
res_24h_vs_NaiveOrdered <- res_24h_vs_Naive[order(res_24h_vs_Naive$pvalue),]


#summary of some basic tallies using the summary function

# summary(res_24h_vs_NaiveOrdered)
# out of 23499 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2661, 11%
# LFC < 0 (down)     : 2286, 9.7%
# outliers [1]       : 0, 0%
# low counts [2]     : 9035, 38%
# (mean count < 12)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


#How many adjusted p-values were less than 0.1?

sum(res_24h_vs_NaiveOrdered$padj < 0.1, na.rm=TRUE)


## 4947
```

### Exploring and exporting results 

```

## MA-plot
#In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable 
#over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored 
#red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as 
#open triangles pointing either up or down.

plotMA(res_24h_vs_Naive, ylim=c(-2,2))
```
![plotMA](https://jl19.githu.io/BASH_Training_Course_2023/Docs/R_Scripts/R_Plots/Rplot_MA_Plot_res_24h_vs_Naive.jpeg)

```
#It is more useful visualize the MA-plot for the shrunken log2 fold changes, which 
#remove the noise associated with log2 fold changes from low count genes without requiring 
#arbitrary filtering thresholds.

plotMA(resNaive_vs_24h_LFC, ylim=c(-2,2))
```

![MA_plot_LFC](https://jl19.github.io/BASH_Training_Course_2023//Docs/R_Scripts/R_Plots/Rplot_MA_plot_resNaive_vs_24h_LFC.jpeg)
```



#### Plot counts    

#It can also be useful to examine the counts of reads for a single gene across the 
#groups. A simple function for making this plot is plotCounts, which normalizes c
#ounts by the estimated size factors (or normalization factors if these were used) and adds a pseudocount of 1/2 to allow for log scale plotting. The counts are grouped by the variables in intgroup, where more than one variable can be specified. Here we specify the gene which had the smallest p value from the results table created above. You can select the gene to plot by rowname or by numeric index.


plotCounts(dds_txi, gene=which.min(res_24h_vs_Naive$padj), intgroup="condition")


#customized plotting


# Save plotcounts to a data frame object
d <- plotCounts(dds_txi, gene=which.min(res_24h_vs_Naive$padj), intgroup="condition", returnData=TRUE)

# Plotting the normalized counts, using the samplenames (rownames of d as labels)
ggplot(d, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  #geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("Comparision of 24h, 2h and Naive") +
  theme(plot.title = element_text(hjust = 0.5))
```
![24h_2h_Navie](https://jl19.github.io/BASH_Training_Course_2023//Docs/R_Scripts/R_Plots/Rplot_Counts_24h_2h_Naive.jpeg)

### Exporting results to CSV files
```
write.csv(as.data.frame(res_24h_vs_NaiveOrdered), 
          file="res_24h_vs_NaiveOrdered_results.csv")


resSig_24h_vs_NaiveOrdered <- subset(res_24h_vs_NaiveOrdered, padj < 0.05)
resSig_24h_vs_NaiveOrdered
```
### Data transformations and visualization
```
#Count data transformations
#In order to test for differential expression, we operate on raw counts and use discrete 
#distributions as described in the previous section on differential expression. However for other 
#downstream analyses – e.g. for visualization or clustering – it might be useful to work with 
#transformed versions of the count data.

#the most obvious choice of transformation is the logarithm. Since count values for a gene can be zero in some conditions (and non-zero in others), some advocate the use of pseudocounts, i.e. transformations of the form:
  
#  y=log2(n+n0)

#where n represents the count values and n0 is a positive constant.


#Extracting transformed values
#These transformation functions return an object of class DESeqTransform which is a subclass of 
#RangedSummarizedExperiment. For ~20 samples, running on a newly created DESeqDataSet, rlog may take 30 seconds, while vst takes less than 1 second. The running times are shorter when using blind=FALSE and if the function DESeq has already been run, because then it is not necessary to re-estimate the dispersion values. The assay function is used to extract the matrix of normalized values.

vsd <- vst(dds_txi, blind=FALSE)
rld <- rlog(dds_txi, blind=FALSE)
head(assay(vsd), 3)

#Effects of transformations on the variance

#the vertical axis in such plots is the square root of the variance over all samples, 
#so including the variance due to the experimental conditions. While a flat curve of the square 
#root of variance over the mean may seem like the goal of such transformations, this may be 
#unreasonable in the case of datasets with many true differences due to the experimental conditions.
# this gives log2(n + 1)

ntd <- normTransform(dds_txi)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

```
### Data quality assessment by sample clustering and visualization
```
#Heatmap of the count matrix
#To explore a count matrix, it is often instructive to look at it as a heatmap. 
#Below we show how to produce such a heatmap for various transformations of the data.

library("pheatmap")
select <- order(rowMeans(counts(dds_txi,normalized=TRUE)),
                decreasing=TRUE)[1:20]
colData(dds_txi)

# the rows/columns of the distance matrix.
# DataFrame with 6 rows and 31 columns
# Run         Age  Assay.Type AvgSpotLen      Bases  BioProject    BioSample
# <character> <character> <character>  <integer>  <numeric> <character>  <character>
#   SRR7457551  SRR7457551    14 weeks     RNA-Seq         75 2244859888 PRJNA450151 SAMN08949698
# SRR7457552  SRR7457552    14 weeks     RNA-Seq         75 2018541831 PRJNA450151 SAMN08949699
# SRR7457555  SRR7457555    14 weeks     RNA-Seq         75 2011980926 PRJNA450151 SAMN08949694
# SRR7457556  SRR7457556    14 weeks     RNA-Seq         75 2029890687 PRJNA450151 SAMN08949695
# SRR7457559  SRR7457559    14 weeks     RNA-Seq         75 1703298074 PRJNA450151 SAMN08949690
# SRR7457560  SRR7457560    14 weeks     RNA-Seq         75 1511125335 PRJNA450151 SAMN08949691
# BioSampleModel     Bytes           Cell_type            Center.Name     Consent
# <character> <integer>         <character>            <character> <character>
#   SRR7457551 Model organism or an.. 883354260 Alveolar Macrophage NORTHWESTERN UNIVERS..      public
# SRR7457552 Model organism or an.. 793650717 Alveolar Macrophage NORTHWESTERN UNIVERS..      public
# SRR7457555 Model organism or an.. 788409431 Alveolar Macrophage NORTHWESTERN UNIVERS..      public
# SRR7457556 Model organism or an.. 796029883 Alveolar Macrophage NORTHWESTERN UNIVERS..      public
# SRR7457559 Model organism or an.. 666357446 Alveolar Macrophage NORTHWESTERN UNIVERS..      public
# SRR7457560 Model organism or an.. 591158034 Alveolar Macrophage NORTHWESTERN UNIVERS..      public
# DATASTORE.filetype DATASTORE.provider   DATASTORE.region  Experiment  Instrument
# <character>        <character>        <character> <character> <character>
#   SRR7457551          sra,fastq              s3,gs s3.us-east-1,gs.US  SRX4328058 NextSeq 500
# SRR7457552          sra,fastq              gs,s3 s3.us-east-1,gs.US  SRX4328057 NextSeq 500
# SRR7457555          fastq,sra              gs,s3 s3.us-east-1,gs.US  SRX4328054 NextSeq 500
# SRR7457556          sra,fastq              s3,gs gs.US,s3.us-east-1  SRX4328053 NextSeq 500
# SRR7457559          sra,fastq              gs,s3 s3.us-east-1,gs.US  SRX4328050 NextSeq 500
# SRR7457560          sra,fastq              s3,gs gs.US,s3.us-east-1  SRX4328049 NextSeq 500
# Isolate   Library.Name LibraryLayout LibrarySelection  LibrarySource     Organism
# <character>    <character>   <character>      <character>    <character>  <character>
#   SRR7457551 24h Reperfusion AM B.. R05_AM_Allo24h        SINGLE              PCR TRANSCRIPTOMIC Mus musculus
# SRR7457552 24h Reperfusion AM B.. R06_AM_Allo24h        SINGLE              PCR TRANSCRIPTOMIC Mus musculus
# SRR7457555 2h Reperfusion AM Bi..  R01_AM_Allo2h        SINGLE              PCR TRANSCRIPTOMIC Mus musculus
# SRR7457556 2h Reperfusion AM Bi..  R02_AM_Allo2h        SINGLE              PCR TRANSCRIPTOMIC Mus musculus
# SRR7457559 Naive AM Biological ..   N01_AM_Naive        SINGLE              PCR TRANSCRIPTOMIC Mus musculus
# SRR7457560 Naive AM Biological ..   N02_AM_Naive        SINGLE              PCR TRANSCRIPTOMIC Mus musculus
# Platform          ReleaseDate    Sample.Name         sex   SRA.Study        strain      tissue
# <character>          <character>    <character> <character> <character>   <character> <character>
#   SRR7457551    ILLUMINA 2018-06-29T00:00:00Z R05_AM_Allo24h        male   SRP151689 Cx3cr1gfp/+B6        Lung
# SRR7457552    ILLUMINA 2018-06-29T00:00:00Z R06_AM_Allo24h        male   SRP151689 Cx3cr1gfp/+B6        Lung
# SRR7457555    ILLUMINA 2018-06-29T00:00:00Z  R01_AM_Allo2h        male   SRP151689 Cx3cr1gfp/+B6        Lung
# SRR7457556    ILLUMINA 2018-06-29T00:00:00Z  R02_AM_Allo2h        male   SRP151689 Cx3cr1gfp/+B6        Lung
# SRR7457559    ILLUMINA 2018-06-29T00:00:00Z   N01_AM_Naive        male   SRP151689 Cx3cr1gfp/+B6        Lung
# SRR7457560    ILLUMINA 2018-06-29T00:00:00Z   N02_AM_Naive        male   SRP151689 Cx3cr1gfp/+B6        Lung
# condition
# <factor>
#   SRR7457551     24h  
# SRR7457552     24h  
# SRR7457555     2h   
# SRR7457556     2h   
# SRR7457559     Naive
# SRR7457560     Naive

#Just select condtion and tissue to from the dds_txi as data.frame
df <- as.data.frame(colData(dds_txi)[,c("condition","tissue")])

#             condition tissue
# SRR7457551       24h   Lung
# SRR7457552       24h   Lung
# SRR7457555        2h   Lung
# SRR7457556        2h   Lung
# SRR7457559     Naive   Lung
# SRR7457560     Naive   Lung

colnames(dds_txi) 
#[1] "SRR7457551" "SRR7457552" "SRR7457555" "SRR7457556" "SRR7457559" "SRR7457560"
#pheatmap A function to draw clustered heatmaps.
#set rownames = colnames
rownames(df) <- colnames(dds_txi)


pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
![Heatmap_Condition_ntd](https://jl19.github.io/BASH_Training_Course_2023/Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_condition_ntd.jpeg)

```
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

```
![Heatmap_Condition_vsd](https://jl19.github.io/BASH_Training_Course_2023/Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_condition_vsd.jpeg)


```
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

![Heatmap_Condition_rld](https://jl19.github.io/BASH_Training_Course_2023//Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_condition_rld.jpeg)

```


#Heatmap of the sample-to-sample distances
#Another use of the transformed data is sample clustering. Here, we apply the dist function
#to the transpose of the transformed count matrix to get sample-to-sample distances.

sampleDists <- dist(t(assay(vsd)))

#A heatmap of this distance matrix gives us an overview over similarities and dissimilarities 
#between samples. We have to provide a hierarchical clustering hc to the heatmap function based on 
#the sample distances, or else the heatmap function would calculate a clustering based on the 
#distances between 
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
rownames(sampleDistMatrix) <- paste(vsd$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

```

![Samples_Distance_Matrix](https://jl19.github.io/BASH_Training_Course_2023//Docs/R_Scripts/R_Plots/Rplot_distance_matrix.jpeg)

### Principal component plot of the samples
```
#Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

#plotPCA(vsd, intgroup=c("condition", "type"))
plotPCA(vsd, intgroup=c("condition"))


pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

```
![PCA_PC1_2](https://jl19.github.io/BASH_Training_Course_2023//Docs/R_Scripts/R_Plots/Rplot_PCA_vsd_condtion_PC1_PC2.jpeg)



### Box plot of the samples 

```
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_txi)[["cooks"]]), range=0, las=2)

```
![Samples_Boxplot](https://jl19.github.io/BASH_Training_Course_2023//Docs/R_Scripts/Results/Rplot_boxplot.jpeg)

```
```

### Tests of log2 fold change above or below a threshold
```
# It is also possible to provide thresholds for constructing Wald tests of significance. Two arguments to the results 
#function allow for threshold-based Wald tests: lfcThreshold, which takes a numeric of a non-negative threshold value, 
#and altHypothesis, which specifies the kind of test. Note that the alternative hypothesis is specified by the user, 
#i.e. those genes which the user is interested in finding, and the test provides p values for the null hypothesis, 
#the complement of the set defined by the alternative. The altHypothesis argument can take one of the following four values, 
#where β is the log2 fold change specified by the name argument, and x is the lfcThreshold.

# 
# greaterAbs - |β|>x - tests are two-tailed
# lessAbs - |β|<x - p values are the maximum of the upper and lower tests
# greater - β>x
# less - β<−x
# The four possible values of altHypothesis are demonstrated in the following code and visually by MA-plots in the following figures.
# 



par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-2.5,2.5)
resGA <- results(dds_txi, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds_txi, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds_txi, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds_txi, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()


```
![Samples_Boxplot](https://jl19.github.io/BASH_Training_Course_2023//Docs/R_Scripts/Results/Rplot_Test_Log2_fc_threshold.jpeg)

[Go to Day2 Practical 6](analyzing-RNA-seq-data-with-DESeq2-data-visualization.md)
## [Go to Day2 Practicals](rna-seq-wes-data-analysis-day2.md/#quantification)


