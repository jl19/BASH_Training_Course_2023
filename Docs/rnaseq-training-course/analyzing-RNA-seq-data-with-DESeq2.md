# Analyzing RNA-seq data with DESeq2 using R

## Data Input and Analysis in R

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
library("vsn")
```
### Set directory
```
dir <- "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583"

setwd("/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/DE_Analysis/")

```
### Read sample information

```
samples <- read.csv("sample_metadata_12samples_stringtie.csv", header = TRUE)
```
### Sample run information
```
samples$Run

#[1] "SRR7457557_control"  "SRR7457558_control"  "SRR7457560_control"  "SRR7457559_control" 
#[5] "SRR7457553_2_hours"  "SRR7457554_2_hours"  "SRR7457555_2_hours"  "SRR7457556_2_hours" 
#[9] "SRR7457561_24_hours" "SRR7457551_24_hours" "SRR7457552_24_hours" "SRR7457562_24_hours"

```
### Set the factor as 24h, 2h and Control replicates 4
```
samples$condition <- factor(c("Control","Control","Control","Control","2h","2h","2h","2h","24h","24h","24h","24h"))
samples$condition
#[1] [1] Control Control Control Control 2h      2h      2h      2h      24h     24h     24h     24h     
#Levels: 24h 2h Control
```


### Loading data from previous step tximport output

Construct a DESeqDataSEt from the txi_stringtie object created in txtimport.R and sample information in samples.
txi is matrices which contain the transcript-level abundance, estimated counts and transcript lengthes etc.

DESeq2 require the dataset has replicates
```
dds_txi_stringtie <- DESeqDataSetFromTximport(txi_stringtie,
                              colData = samples,
                              design = ~ condition)
```
## Differential expression analysis

```
dds_txi_stringtie <- DESeq(dds_txi_stringtie)

# estimating size factors
# using 'avgTxLength' from assays(dds), correcting for library size
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
```

### Result for Control vs 24h

```
res_stringtie <- results(dds_txi_stringtie)

res_stringtie

# log2 fold change (MLE): condition Control vs 24h 
# Wald test p-value: condition Control vs 24h 
# DataFrame with 48321 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#        baseMean log2FoldChange     lfcSE        stat      pvalue        padj
#<numeric>      <numeric> <numeric>   <numeric>   <numeric>   <numeric>
# 0610006L08Rik   0.00000             NA        NA          NA          NA          NA
#0610007P14Rik 373.64780     -1.4565690  0.133412 -10.9178450 9.47165e-28 4.61882e-26
#0610009B22Rik 186.75776     -0.0907768  0.120443  -0.7536938 4.51033e-01 6.14060e-01
#0610009E02Rik   1.52217     -0.1616363  2.198465  -0.0735224 9.41390e-01          NA
#0610009L18Rik  16.35806      0.8412978  0.370995   2.2676778 2.33489e-02 6.16440e-02
#...                 ...            ...       ...         ...         ...         ...
#Zyg11a             0.00             NA        NA          NA          NA          NA
#Zyg11b          1084.92      0.2755864 0.0645755    4.267661 1.97533e-05 1.13754e-04
#Zyx             4717.40      0.9214825 0.0926328    9.947691 2.58102e-23 9.76610e-22
#Zzef1           2310.02      0.0700811 0.0785500    0.892184 3.72295e-01 5.38625e-01
#Zzz3            1141.16      0.3587183 0.1220460    2.939205 3.29055e-03 1.14792e-02

```
Note that we could have specified the coefficient or contrast we want to build a results table for, using either of the following equivalent commands:

### Check existing names of the object dds_txi_stringtie
```
resultsNames(dds_txi_stringtie)

#[1] "Intercept"              "condition_2h_vs_24h"    "condition_Control_vs_24h"
```
### Compare two conditions using name
```
res_control_vs_24h <- results(dds_txi_stringtie, name=c("condition_Control_vs_24h")
res_control_vs_24h

#log2 fold change (MLE): condition Control vs 24h 
#Wald test p-value: condition Control vs 24h 
#DataFrame with 48321 rows and 6 columns
#baseMean log2FoldChange     lfcSE        stat      pvalue        padj
#<numeric>      <numeric> <numeric>   <numeric>   <numeric>   <numeric>
# 0610006L08Rik   0.00000             NA        NA          NA          NA          NA
#0610007P14Rik 373.64780     -1.4565690  0.133412 -10.9178450 9.47165e-28 4.61882e-26
#0610009B22Rik 186.75776     -0.0907768  0.120443  -0.7536938 4.51033e-01 6.14060e-01
#0610009E02Rik   1.52217     -0.1616363  2.198465  -0.0735224 9.41390e-01          NA
#0610009L18Rik  16.35806      0.8412978  0.370995   2.2676778 2.33489e-02 6.16440e-02
#...                 ...            ...       ...         ...         ...         ...
#Zyg11a             0.00             NA        NA          NA          NA          NA
#Zyg11b          1084.92      0.2755864 0.0645755    4.267661 1.97533e-05 1.13754e-04
```
### Compare two conditions using contrast
Contrast is the factor set in txiimport.R which are Control, 2h and 24h

```
res_control_vs_2h <- results(dds_txi_stringtie, contrast=c("condition","Control","2h")
res_24h_vs_2h <- results(dds_txi_stringtie, contrast = c("condition", "24h", "2h"))
res_24h_vs_Control <- results(dds_txi_stringtie, contrast=c("condition","24h", "Control"))
res_2h_vs_Control <- results(dds_txi_stringtie, contrast=c("condition","2h", "Control"))
```

### Log fold change shrinkage for visualization and ranking


Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. 
To shrink the LFC, we pass the dds object to the function lfcShrink. Below we specify to use the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.

DESeq2 provides the dds object and the name or number of the coefficient we want to shrink, where the number refers to the order of the coefficient as it appears in resultsNames(dds_txi_stringtie)

### List resultsNames
```
resultsNames(dds_txi_stringtie)
## [1] "Intercept"                     "condition_2h_vs_24h"    "condition_Control_vs_24h"

resControl_vs_24h_LFC <- lfcShrink(dds_txi_stringtie, coef="condition_Control_vs_24h", type="apeglm")
```
Using 'apeglm' for LFC shrinkage. If used in published research, please cite:Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895
```
resControl_vs_24h_LFC

# log2 fold change (MAP): condition Control vs 24h 
#Wald test p-value: condition Control vs 24h 
#DataFrame with 48321 rows and 5 columns
# baseMean log2FoldChange     lfcSE      pvalue        padj
# <numeric>      <numeric> <numeric>   <numeric>   <numeric>
#   0610006L08Rik   0.00000             NA        NA          NA          NA
# 0610007P14Rik 373.64780    -1.43484405  0.134099 9.47165e-28 4.61882e-26
# 0610009B22Rik 186.75776    -0.08307114  0.115533 4.51033e-01 6.14060e-01
# 0610009E02Rik   1.52217    -0.00747285  0.387772 9.41390e-01          NA
# 0610009L18Rik  16.35806     0.59583003  0.375758 2.33489e-02 6.16440e-02
# ...                 ...            ...       ...         ...         ...
# Zyg11a             0.00             NA        NA          NA          NA
# Zyg11b          1084.92      0.2698079 0.0641510 1.97533e-05 1.13754e-04
# Zyx             4717.40      0.9078606 0.0929471 2.58102e-23 9.76610e-22
# Zzef1           2310.02      0.0687265 0.0771193 3.72295e-01 5.38625e-01
# Zzz3            1141.16      0.3352972 0.1200833 3.29055e-03 1.14792e-02

```

### Ordering results using p-values and adjusted p-values 
```
res_24h_vs_Control_Ordered <- res_24h_vs_Control[order(res_24h_vs_Control$pvalue),]
res_24h_vs_Control_Ordered

# log2 fold change (MLE): condition 24h vs Control 
# Wald test p-value: condition 24h vs Control 
# DataFrame with 48321 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#   2810417H13Rik   1368.76        2.95878 0.0898318   32.9369 6.52038e-238 1.18919e-233
# Pcna            4360.31        1.46259 0.0471930   30.9918 6.96217e-211 6.34880e-207
# Rrm1            2486.12        2.17967 0.0737897   29.5390 9.10047e-192 5.53248e-188
# Rrm2            2041.04        2.58360 0.0891632   28.9760 1.31903e-184 6.01413e-181
# Nusap1           984.91        2.29597 0.0834374   27.5173 1.09074e-166 3.97860e-163
# ...                 ...            ...       ...       ...          ...          ...
# Zscan4e               0             NA        NA        NA           NA           NA
# Zscan4f               0             NA        NA        NA           NA           NA
# Zscan5b               0             NA        NA        NA           NA           NA
```
### Summary of some basic tallies using the summary function
```
summary(res_24h_vs_Control_Ordered)
# out of 25261 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3873, 15%
# LFC < 0 (down)     : 3771, 15%
# outliers [1]       : 30, 0.12%
# low counts [2]     : 6993, 28%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
```

### How many adjusted p-values were less than 0.05?
```
sum(res_24h_vs_Control_Ordered$padj < 0.01, na.rm=TRUE)

#[1]5128
```
## Data transformations and visualization      

### MA-plot
In DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples in the DESeqDataSet. Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down.

```
plotMA(res_24h_vs_Naive, ylim=c(-2,2))
```
![plotMA](https://jl19.github.io/BASH_Training_Course_2023/Docs/R_Scripts/R_Plots/Rplot_MA_Plot_res_24h_vs_Naive.jpeg)


It is more useful visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
```
plotMA(resNaive_vs_24h_LFC, ylim=c(-2,2))


![MA_plot_LFC](https://jl19.github.io/BASH_Training_Course_2023//Docs/R_Scripts/R_Plots/Rplot_MA_plot_resNaive_vs_24h_LFC.jpeg)
```

### Plot Counts    

#Plot Counts

It can also be useful to examine the counts of reads for a single gene across the groups. A simple function for making this plot is plotCounts, which normalizes counts by the estimated size factors (or normalization factors if these were used) and adds a pseudocount of 1/2 to allow for log scale plotting. The counts are grouped by the variables in intgroup, where more than one variable can be specified. Here we specify the gene which had the smallest p value from the results table created above. You can select the gene to plot by rowname or by numeric index.

```
#The output below gives out the list of transcript/gene sorted according to adjusted p value

# log2 fold change (MLE): condition 24h vs Control 
# Wald test p-value: condition 24h vs Control 
# DataFrame with 48321 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
# 2810417H13Rik   1368.76        2.95878 0.0898318   32.9369 6.52038e-238 1.18919e-233
# Pcna            4360.31        1.46259 0.0471930   30.9918 6.96217e-211 6.34880e-207
# Rrm1            2486.12        2.17967 0.0737897   29.5390 9.10047e-192 5.53248e-188
# Rrm2            2041.04        2.58360 0.0891632   28.9760 1.31903e-184 6.01413e-181
# Nusap1           984.91        2.29597 0.0834374   27.5173 1.09074e-166 3.97860e-163

# "which.min" gives gene with the lowest p-value

plotCounts(dds_txi_stringtie, gene=which.min(res_24h_vs_Control$padj), intgroup="condition")

#Use gene ID "2810417H13Rik" to plot the counts
plotCounts(dds_txi_stringtie, "2810417H13Rik", intgroup="condition")
#Use gene name Pcna to plot the counts
plotCounts(dds_txi_stringtie, "Pcna", intgroup="condition")

#customized plotting

#Save plotCounts to a data frame object

df_Pcna <- plotCounts(dds_txi_stringtie,"Pcna" , intgroup="condition", returnData=TRUE)


# Plotting the normalized counts
ggplot(df_Pcna, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  theme_bw() +
  ggtitle("Proliferating cell nuclear antigen (PCNA)") +
  theme(plot.title = element_text(hjust = 0.5))
```
![24h_2h_Navie](https://jl19.github.io/BASH_Training_Course_2023//Docs/R_Scripts/R_Plots/Rplot_Counts_24h_2h_Naive.jpeg)

### Exporting results to CSV files
```
write.csv(as.data.frame(res_24h_vs_Control_Ordered), 
          file="res_24h_vs_Control_Ordered_results.csv")

#Extract subset of the dataset
resSig_24h_vs_Control_Ordered <- subset(res_24h_vs_Control_Ordered, padj < 0.05)
resSig_24h_vs_Control_Ordered

# log2 fold change (MLE): condition 24h vs Control 
# Wald test p-value: condition 24h vs Control 
# DataFrame with 6666 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#   2810417H13Rik    1368.76        2.95878 0.0898318   32.9369 6.52038e-238 1.18919e-233
# Pcna             4360.31        1.46259 0.0471930   30.9918 6.96217e-211 6.34880e-207
# Rrm1             2486.12        2.17967 0.0737897   29.5390 9.10047e-192 5.53248e-188
# Rrm2             2041.04        2.58360 0.0891632   28.9760 1.31903e-184 6.01413e-181
# Nusap1            984.91        2.29597 0.0834374   27.5173 1.09074e-166 3.97860e-163

```
## Data Transformations and Visualization

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

```
vsd <- vst(dds_txi, blind=FALSE)
rld <- rlog(dds_txi, blind=FALSE)
head(assay(vsd), 3)
```
#Effects of transformations on the variance

#the vertical axis in such plots is the square root of the variance over all samples, 
#so including the variance due to the experimental conditions. While a flat curve of the square 
#root of variance over the mean may seem like the goal of such transformations, this may be 
#unreasonable in the case of datasets with many true differences due to the experimental conditions.
# this gives log2(n + 1)
```
ntd <- normTransform(dds_txi)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

```
### Data quality assessment by sample clustering and visualization

#Heatmap of the count matrix
#To explore a count matrix, it is often instructive to look at it as a heatmap. 
#Below we show how to produce such a heatmap for various transformations of the data.
```
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
```
Just select condtion and tissue to from the dds_txi as data.frame

```
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


