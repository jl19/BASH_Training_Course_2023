# Data transformations and visualization

Continued from previous DE analysis

## Count data transformations
```
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
![Heatmap_Condition_ntd](/Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_condition_ntd.jpeg)

```
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```
![Heatmap_Condition_vsd](/Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_condition_vsd.jpeg)
```
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
```

![Heatmap_Condition_rld](/Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_condition_rld.jpeg)
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
![Samples_Distance_Matrix](/Docs/R_Scripts/R_Plots/Rplot_distance_matrix.jpeg)

#### Principal component plot of the samples

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
![PCA_PC1_2](/Docs/R_Scripts/R_Plots/Rplot_PCA_vsd_condtion_PC1_PC2.jpeg)

#### Box plot of the samples

```
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_txi)[["cooks"]]), range=0, las=2)

```
![Samples_Boxplot](/R_Scripts/Results/Rplot_boxplot.jpeg)

```
```
#### Tests of log2 fold change above or below a threshold

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
![Samples_Boxplot](/Docs/R_Scripts/Results/Rplot_Test_Log2_fc_threshold.jpeg)


## [Go to Day2 Practicals](rna-seq-wes-data-analysis-day2/#quantification)


