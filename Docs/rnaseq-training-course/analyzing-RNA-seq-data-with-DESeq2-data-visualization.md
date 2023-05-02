## Data Transformations and Visualization

### Count data transformations

In order to test for differential expression, raw counts and use discrete distributions as described in the previous section on differential expression. However for other downstream analyses – e.g. for visualization or clustering – it might be useful to work with transformed versions of the count data.

The most obvious choice of transformation is the logarithm. Since count values for a gene can be zero in some conditions (and non-zero in others), some advocate the use of pseudocounts, i.e. transformations of the form:
  
y=log2(n+n0)

where n represents the count values and n0 is a positive constant.


### Extracting transformed values
These transformation functions return an object of class DESeqTransform which is a subclass of RangedSummarizedExperiment. For ~20 samples, running on a newly created DESeqDataSet, rlog may take 30 seconds, while vst takes less than 1 second. The running times are shorter when using blind=FALSE and if the function DESeq has already been run, because then it is not necessary to re-estimate the dispersion values. The assay function is used to extract the matrix of normalized values.

```
vsd <- vst(dds_txi, blind=FALSE)
rld <- rlog(dds_txi, blind=FALSE)
head(assay(vsd), 3)
```
### Effects of transformations on the variance

The vertical axis in such plots is the square root of the variance over all samples, so including the variance due to the experimental conditions. While a flat curve of the square root of variance over the mean may seem like the goal of such transformations, this may be unreasonable in the case of datasets with many true differences due to the experimental conditions. this gives log2(n + 1)
```
ntd <- normTransform(dds_txi)

meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))

```
### Data quality assessment by sample clustering and visualization

#### Heatmap of the count matrix
Heatmap is useful way to explore a count matrix

```
select <- order(rowMeans(counts(dds_txi_stringtie,normalized=TRUE)),
                decreasing=TRUE)[1:20]


colData(dds_txi_stringtie)

### the rows/columns of the distance matrix.


DataFrame with 12 rows and 33 columns
                                    Run         Age  Assay.Type AvgSpotLen      Bases
                            <character> <character> <character>  <integer>  <numeric>
SRR7457557_control   SRR7457557_control    14 weeks     RNA-Seq         75 1441343365
SRR7457558_control   SRR7457558_control    14 weeks     RNA-Seq         75 1178449406
SRR7457559_control   SRR7457559_control    14 weeks     RNA-Seq         75 1703298074
SRR7457560_control   SRR7457560_control    14 weeks     RNA-Seq         75 1511125335
SRR7457553_2_hours   SRR7457553_2_hours    14 weeks     RNA-Seq         75 1705376082
...                                 ...         ...         ...        ...        ...
SRR7457556_2_hours   SRR7457556_2_hours    14 weeks     RNA-Seq         75 2029890687
SRR7457551_24_hours SRR7457551_24_hours    14 weeks     RNA-Seq         75 2244859888
SRR7457552_24_hours SRR7457552_24_hours    14 weeks     RNA-Seq         75 2018541831
SRR7457562_24_hours SRR7457561_24_hours    14 weeks     RNA-Seq         75 1786917960
SRR7457561_24_hours SRR7457562_24_hours    14 weeks     RNA-Seq         75 1730344657
```
### Just select condition and tissue to from the dds_txi_stringtie as data.frame

```
df <- as.data.frame(colData(dds_txi_stringtie)[,c("condition","TISSUE")])

# condition TISSUE
# SRR7457557_control    Control   Lung
# SRR7457558_control    Control   Lung
# SRR7457559_control    Control   Lung
# SRR7457560_control    Control   Lung
# SRR7457553_2_hours         2h   Lung
# SRR7457554_2_hours         2h   Lung
# SRR7457555_2_hours         2h   Lung
# SRR7457556_2_hours         2h   Lung
# SRR7457551_24_hours       24h   Lung
# SRR7457552_24_hours       24h   Lung
# SRR7457562_24_hours       24h   Lung
# SRR7457561_24_hours       24h   Lung
```
```
colnames(dds_txi_stringtie) 
# [1] "SRR7457551_24_hours" "SRR7457552_24_hours" "SRR7457553_2_hours" 
#[4] "SRR7457554_2_hours"  "SRR7457555_2_hours"  "SRR7457556_2_hours" 
#[7] "SRR7457557_control"  "SRR7457558_control"  "SRR7457560_control" 
#[10] "SRR7457562_24_hours" "SRR7457559_control"  "SRR7457561_24_hours"
```

### Pheatmap A function to draw clustered heatmaps.
```
#set rownames(df_stringtie) = colnames(dds_txi_stringtie)
rownames(df_stringtie) <- colnames(dds_txi_stringtie)


pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_stringtie)
```
![Heatmap_Condition_ntd](https://jl19.github.io/BASH_Training_Course_2023/Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_condition_ntd.png)
```
```
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_stringtie)

```
```
![Heatmap_Condition_vsd](https://jl19.github.io/BASH_Training_Course_2023/Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_condition_vsd.png)


```
```
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_stringtie)
```
```
![Heatmap_Condition_rld](https://jl19.github.io/BASH_Training_Course_2023/Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_condition_rld.png)

```
```
### Heatmap of the sample-to-sample distances
Another use of the transformed data is sample clustering. Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.

#### Obtain the sample euclidean distances

```
sampleDists <- dist(t(assay(vsd)))
```
A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples. We have to provide a hierarchical clustering hc to the heatmap function based on the sample distances, or else the heatmap function would calculate a clustering based on the distances between 
```
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
```
#### Add names based on condition
```
rownames(sampleDistMatrix) <- paste(vsd$condition)
colnames(sampleDistMatrix) <- NULL
## Define colors to use for the heatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

## Make the heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```
![Heatmap_sample_distribution_vsd](https://jl19.github.io/BASH_Training_Course_2023/Docs/R_Scripts/R_Plots/Rplot_heatmap_sample_distribution_vsd.png)

This plot shows how samples are clustered based on their euclidean distance using the regularized log transformed count data. This figure gives an overview of how the samples are hierarchically clustered. It is a complementary figure to the PCA plot.


### Principal component plot of the samples

Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

```
plotPCA(vsd, intgroup=c("condition"))

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

```


## [Go to Day2 Practicals](rna-seq-wes-data-analysis-day2.md/#quantification)


