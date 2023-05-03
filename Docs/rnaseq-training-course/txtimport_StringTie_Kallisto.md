# Import Transcript-level Estimates

> - Imports transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages. 
> - Average transcript length, weighted by sample-specific transcript abundance estimates, is provided as a matrix which can be used as an offset for different expression of gene-level counts.


## Import from Kallisto Output

### Install required R packages

``` 
if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install("ballgown")
BiocManager::install("rhdf5")
BiocManager::install("tximport")
BiocManager::install("haven")
``` 
### Install required library
``` 
library(dplyr) # data wrangling
library(ggplot2) # plotting
library(DESeq2) # rna-seq
library(edgeR) # rna-seq
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5) # read/convert kalisto output files.  
``` 

### Set working directory
``` 
dir <- "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583"

setwd("/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/DE_Analysis/")

``` 
### Read sample information
``` 
samples <- read.csv("sample_metadata_12samples_kallisto.csv", header = TRUE)
``` 
### Sample run information

``` 
samples
                           Run      Age Assay.Type AvgSpotLen
1   SRR7457557_control-trimmed 14 weeks    RNA-Seq         75
2   SRR7457558_control-trimmed 14 weeks    RNA-Seq         75
3   SRR7457559_control-trimmed 14 weeks    RNA-Seq         75
4   SRR7457560_control-trimmed 14 weeks    RNA-Seq         75
5   SRR7457553_2_hours-trimmed 14 weeks    RNA-Seq         75
6   SRR7457554_2_hours-trimmed 14 weeks    RNA-Seq         75
7   SRR7457555_2_hours-trimmed 14 weeks    RNA-Seq         75
8   SRR7457556_2_hours-trimmed 14 weeks    RNA-Seq         75
9  SRR7457551_24_hours-trimmed 14 weeks    RNA-Seq         75
10 SRR7457552_24_hours-trimmed 14 weeks    RNA-Seq         75
11 SRR7457561_24_hours-trimmed 14 weeks    RNA-Seq         75
12 SRR7457562_24_hours-trimmed 14 weeks    RNA-Seq         75
        Bases  BioProject    BioSample           BioSampleModel
1  1441343365 PRJNA450151 SAMN08949692 Model organism or animal
2  1178449406 PRJNA450151 SAMN08949693 Model organism or animal
3  1703298074 PRJNA450151 SAMN08949690 Model organism or animal
4  1511125335 PRJNA450151 SAMN08949691 Model organism or animal
5  1705376082 PRJNA450151 SAMN08949696 Model organism or animal
``` 
To access one variable in a dataset, use the dollar sign “$”. For example, $Run returns the Run variable (the Run column).
``` 
samples$Run

# [1] "SRR7457557_control-trimmed"  "SRR7457558_control-trimmed"  "SRR7457559_control-trimmed" 
# [4] "SRR7457560_control-trimmed"  "SRR7457553_2_hours-trimmed"  "SRR7457554_2_hours-trimmed" 
# [7] "SRR7457555_2_hours-trimmed"  "SRR7457556_2_hours-trimmed"  "SRR7457551_24_hours-trimmed"
# [10] "SRR7457552_24_hours-trimmed" "SRR7457561_24_hours-trimmed" "SRR7457562_24_hours-trimmed"
``` 
### Set the path
``` 
file_path_kallisto <- file.path(dir,"kallisto_output", samples$Run)
file_path_kallisto 

#[1] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457557_control-trimmed" 
# [2] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457558_control-trimmed" 
# [3] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457559_control-trimmed" 
# [4] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457560_control-trimmed" 
# [5] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457553_2_hours-trimmed" 
# [6] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457554_2_hours-trimmed" 
# [7] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457555_2_hours-trimmed" 
# [8] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457556_2_hours-trimmed" 
# [9] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457551_24_hours-trimmed"
#[10] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457552_24_hours-trimmed"
#[11] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457561_24_hours-trimmed"
#[12] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457562_24_hours-trimmed"
``` 
### List files in the path
``` 
list.files(file_path_kallisto)

#[1] "abundance.h5"  "abundance.h5"  "abundance.h5"  "abundance.h5" 
#[5] "abundance.h5"  "abundance.h5"  "abundance.h5"  "abundance.h5" 
#[9] "abundance.h5"  "abundance.h5"  "abundance.h5"  "abundance.h5" 
#[13] "abundance.tsv" "abundance.tsv" "abundance.tsv" "abundance.tsv"
#[17] "abundance.tsv" "abundance.tsv" "abundance.tsv" "abundance.tsv"
#[21] "abundance.tsv" "abundance.tsv" "abundance.tsv" "abundance.tsv"
#[25] "run_info.json" "run_info.json" "run_info.json" "run_info.json"
#[29] "run_info.json" "run_info.json" "run_info.json" "run_info.json"
#[33] "run_info.json" "run_info.json" "run_info.json" "run_info.json"
``` 
### Assign files from kallisto output
``` 
files_kallisto <- file.path(dir, "kallisto_output", samples$Run, "abundance.tsv")
files_kallisto
# [1] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457557_control-trimmed/abundance.tsv" 
# [2] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457558_control-trimmed/abundance.tsv" 
# [3] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457559_control-trimmed/abundance.tsv" 
# [4] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457560_control-trimmed/abundance.tsv" 
# [5] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457553_2_hours-trimmed/abundance.tsv" 
# [6] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457554_2_hours-trimmed/abundance.tsv" 
# [7] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457555_2_hours-trimmed/abundance.tsv" 
# [8] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457556_2_hours-trimmed/abundance.tsv" 
# [9] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457551_24_hours-trimmed/abundance.tsv"
# [10] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457552_24_hours-trimmed/abundance.tsv"
# [11] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457561_24_hours-trimmed/abundance.tsv"
# [12] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457562_24_hours-trimmed/abundance.tsv"

``` 
### Check if the files exist
``` 
file.exists(files_kallisto)

#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
``` 
### Assign sample names to files
``` 
names(files_kallisto) <- c('SRR7457557_control','SRR7457558_control','SRR7457559_control','SRR7457560_control','SRR7457553_2_hours','SRR7457554_2_hours','SRR7457555_2_hours','SRR7457556_2_hours','SRR7457551_24_hours','SRR7457552_24_hours','SRR7457561_24_hours','SRR7457562_24_hours')

names(files_kallisto)
# [1] "SRR7457557_control"  "SRR7457558_control"  "SRR7457559_control"  "SRR7457560_control"  "SRR7457553_2_hours" 
# [6] "SRR7457554_2_hours"  "SRR7457555_2_hours"  "SRR7457556_2_hours"  "SRR7457551_24_hours" "SRR7457552_24_hours"
# [11] "SRR7457561_24_hours" "SRR7457562_24_hours"
``` 
### List the sample information
``` 
samples

#          Run      Age Assay.Type AvgSpotLen      Bases  BioProject    BioSample
# 1  SRR7457551_24_hours-trimmed 14 weeks    RNA-Seq         75 2244859888 PRJNA450151 SAMN08949698
# 2  SRR7457552_24_hours-trimmed 14 weeks    RNA-Seq         75 2018541831 PRJNA450151 SAMN08949699
# 3   SRR7457553_2_hours-trimmed 14 weeks    RNA-Seq         75 1705376082 PRJNA450151 SAMN08949696
# 4   SRR7457554_2_hours-trimmed 14 weeks    RNA-Seq         75 1414060046 PRJNA450151 SAMN08949697
# 5   SRR7457555_2_hours-trimmed 14 weeks    RNA-Seq         75 2011980926 PRJNA450151 SAMN08949694
# 6   SRR7457556_2_hours-trimmed 14 weeks    RNA-Seq         75 2029890687 PRJNA450151 SAMN08949695
# 7   SRR7457557_control-trimmed 14 weeks    RNA-Seq         75 1441343365 PRJNA450151 SAMN08949692
# 8   SRR7457558_control-trimmed 14 weeks    RNA-Seq         75 1178449406 PRJNA450151 SAMN08949693
# 9   SRR7457560_control-trimmed 14 weeks    RNA-Seq         75 1511125335 PRJNA450151 SAMN08949691
# 10 SRR7457562_24_hours-trimmed 14 weeks    RNA-Seq         75 1730344657 PRJNA450151 SAMN08949701
# 11  SRR7457559_control-trimmed 14 weeks    RNA-Seq         75 1703298074 PRJNA450151 SAMN08949690
# 12 SRR7457561_24_hours-trimmed 14 weeks    RNA-Seq         75 1786917960 PRJNA450151 SAMN08949700
```      
### List the run names
``` 
samples$Run

#[1] "SRR7457551_24_hours-trimmed" "SRR7457552_24_hours-trimmed" "SRR7457553_2_hours-trimmed" 
#[4] "SRR7457554_2_hours-trimmed"  "SRR7457555_2_hours-trimmed"  "SRR7457556_2_hours-trimmed" 
#[7] "SRR7457557_control-trimmed"  "SRR7457558_control-trimmed"  "SRR7457560_control-trimmed" 
#[10] "SRR7457562_24_hours-trimmed" "SRR7457559_control-trimmed"  "SRR7457561_24_hours-trimmed"

``` 
## Creating tx2gene dataframe

Transcripts need to be associated with gene IDs for gene-level summarization. If that information is present in the files, we can skip this step. For Salmon, Sailfish, and kallisto the files only provide the transcript ID. We first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. The column names do not matter but this column order must be used. The transcript ID must be the same one used in the abundance files.

Creating this tx2gene data.frame can be accomplished from a TxDb object and the select function from the AnnotationDbi package. The following code could be used to construct such a table:
  
For kallisto the files only provide the transcript ID. We first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. The column names do not matter but this column order must be used. The transcript ID must be the same one used in the abundance files.

Creating this tx2gene data.frame can be accomplished from a TxDb object and the select function from the AnnotationDbi package. The following code could be used to construct such a table:

### Load TxDb annotation package
``` 
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", force = TRUE)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#It is R interface to the databases

##load thel ibrary
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 

##list the contents that areloaded intomemory 
ls('package:TxDb.Mmusculus.UCSC.mm10.knownGene') 

##show the db object that is loaded by calling it's name 
TxDb.Mmusculus.UCSC.mm10.knownGene

##assign the db object to txdb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

txdb
#TxDb object:
  # Db type: TxDb
  # Supporting package: GenomicFeatures
  # Data source: UCSC
  # Genome: mm10
  # Organism: Mus musculus
  # Taxonomy ID: 10090
  # UCSC Table: knownGene
  # UCSC Track: GENCODE VM23
  # Resource URL: http://genome.ucsc.edu/
  # Type of Gene ID: Entrez Gene ID
  # Full dataset: yes
# miRBase build ID: NA
# transcript_nrow: 142446
# exon_nrow: 447558
# cds_nrow: 243727
# Db created by: GenomicFeatures package from Bioconductor
# Creation time: 2019-10-21 20:52:26 +0000 (Mon, 21 Oct 2019)
# GenomicFeatures version at creation time: 1.37.4
# RSQLite version at creation time: 2.1.2
# DBSCHEMAVERSION: 1.2

Count the number of genes and transcripts

``` 
genelen <- length(genes(TxDb.Mmusculus.UCSC.mm10.knownGene))

#66 genes were dropped because they have exons located on both strands of the
#same reference sequence or on more than one reference sequence, so cannot be
#represented by a single genomic range.
#Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList
#object, or use suppressMessages() to suppress this message.
#[1] 24528

``` 

columns(txdb)
# [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSPHASE"   "CDSSTART"   "CDSSTRAND"  "EXONCHROM"  "EXONEND"    "EXONID"     "EXONNAME"  
# [12] "EXONRANK"   "EXONSTART"  "EXONSTRAND" "GENEID"     "TXCHROM"    "TXEND"      "TXID"       "TXNAME"     "TXSTART"    "TXSTRAND"   "TXTYPE"    

keytypes(txdb) 

## [1] "CDSID"    "CDSNAME"  "EXONID"   "EXONNAME" "GENEID"   "TXID"     "TXNAME"  

# using ‘exonsBy’ or ‘cdsBy’ with ‘by = "tx"’, the returned exons or CDS are ordered by ascending rank for each transcript,

mm10_transcripts <- exonsBy(txdb, by="tx", use.names=TRUE) 
mm10_transcripts 

#GRangesList object of length 142446:
# $ENSMUST00000193812.1
# GRanges object with 1 range and 3 metadata columns:
#   seqnames          ranges strand |   exon_id   exon_name exon_rank
# <Rle>       <IRanges>  <Rle> | <integer> <character> <integer>
#   [1]     chr1 3073253-3074322      + |         1        <NA>         1
# -------
#   seqinfo: 66 sequences (1 circular) from mm10 genome
# 
# $ENSMUST00000082908.1
# GRanges object with 1 range and 3 metadata columns:
#   seqnames          ranges strand |   exon_id   exon_name exon_rank
# <Rle>       <IRanges>  <Rle> | <integer> <character> <integer>
#   [1]     chr1 3102016-3102125      + |         2        <NA>         1
# -------
#   seqinfo: 66 sequences (1 circular) from mm10 genome
# 
# $ENSMUST00000192857.1
# GRanges object with 1 range and 3 metadata columns:
#   seqnames          ranges strand |   exon_id   exon_name exon_rank
# <Rle>       <IRanges>  <Rle> | <integer> <character> <integer>
#   [1]     chr1 3252757-3253236      + |         3        <NA>         1
# -------
#   seqinfo: 66 sequences (1 circular) from mm10 genome
# 
# ...
# <142443 more elements>

k <- keys(txdb, keytype = "TXNAME")

# [1] "ENSMUST00000193812.1"  "ENSMUST00000082908.1"  "ENSMUST00000192857.1"  "ENSMUST00000161581.1"  "ENSMUST00000192183.1"  "ENSMUST00000193244.1" 
# [7] "ENSMUST00000194454.1"  "ENSMUST00000193450.1"  "ENSMUST00000194935.1"  "ENSMUST00000195361.1"  "ENSMUST00000180019.1"  "ENSMUST00000192738.1" 
# [13] "ENSMUST00000193658.1"  "ENSMUST00000134384.7"  "ENSMUST00000027036.10" "ENSMUST00000150971.7"  "ENSMUST00000155020.1"  "ENSMUST00000119612.8" 
# [19] "ENSMUST00000137887.7"  "ENSMUST00000115529.7"  "ENSMUST00000131119.1"  "ENSMUST00000141278.1"  "ENSMUST00000081551.13" "ENSMUST00000165720.2" 
# [25] "ENSMUST00000144339.1"  "ENSMUST00000169520.1"  "ENSMUST00000192847.5"  "ENSMUST00000044369.12" "ENSMUST00000194676.5"  "ENSMUST00000194301.5" 
# [31] "ENSMUST00000194978.5"  "ENSMUST00000192029.5"  "ENSMUST00000192698.2"  "ENSMUST00000192142.1"  "ENSMUST00000193241.1"  "ENSMUST00000194850.1" 

tx2gene_db <- select(txdb, keys = k, columns="GENEID", keytype = "TXNAME") 

select(txdb, keys = k, columns="GENEID", keytype = "TXNAME") 


# replace NA in the tx2gene_db with blank

tx2gene_db_omitNA <- replace(tx2gene_db, is.na(tx2gene_db), "")
tx2gene_db_omitNA

# TXNAME                GENEID
# 1    ENSMUST00000193812.1  ENSMUSG00000102693.1
# 2    ENSMUST00000082908.1  ENSMUSG00000064842.1
# 3    ENSMUST00000192857.1  ENSMUSG00000102851.1
# 4    ENSMUST00000161581.1  ENSMUSG00000089699.1
# 5    ENSMUST00000192183.1  ENSMUSG00000103147.1
# 6    ENSMUST00000193244.1  ENSMUSG00000102348.1

#Import transcript-level abundance and counts for transcript
#Different type of data generated by different software can be imported ("none", "salmon", "sailfish", "alevin", "kallisto", "rsem", "stringtie")
#########################################################################################################################
#Here is the import type is kallisto
txi_kallisto <- tximport(files_kallisto, type = "kallisto", tx2gene = tx2gene_db_omitNA, ignoreAfterBar = TRUE,dropInfReps = TRUE)

# 1 2 3 4 5 6 
# transcripts missing from tx2gene: 2122
# summarizing abundance
# summarizing counts
# summarizing length

#Summary of the txi
summary(txi_kallisto)

#                     Length Class  Mode     
# abundance           252864  -none- numeric  
# counts              252864  -none- numeric  
# length              252864  -none- numeric  
# countsFromAbundance      1 -none- character

#Take the counts data from the output from tximport
countData_kallisto <- txi_kallisto$counts

countData_kallisto
# SRR7457551_24_hours SRR7457552_24_hours SRR7457553_2_hours SRR7457554_2_hours
#               9.899398e+05        9.752320e+05       8.422768e+05       6.823460e+05
# 100009600        1.903941e+01        2.174964e+01       1.623218e+01       2.352240e+00
# 100009609        1.138890e+01        4.314980e+00       7.382950e+00       7.280110e+00
# 100009614        0.000000e+00        0.000000e+00       0.000000e+00       0.000000e+00
# 100012           0.000000e+00        0.000000e+00       0.000000e+00       0.000000e+00
# 100017           1.495003e+03        1.263000e+03       1.720005e+03       1.088000e+03
# 100019           2.076546e+03        2.019474e+03       1.960450e+03       1.199829e+03
``` 
### Write the counts to a file
``` 
write.csv(countData_kallisto,file="countData_kallisto.csv")
``` 
## Import transcript-level estimates from Stringtie output

Imports transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages. Average transcript length, weighted by sample-specific transcript abundance estimates, is provided as a matrix which can be used as an offset for different expression of gene-level counts.
``` 
if (!requireNamespace("BiocManager", quietly=TaRUE))
install.packages("BiocManager")
BiocManager::install("ballgown")
BiocManager::install("tximport")
BiocManager::install("haven")

library("tximport")
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
### List sample run information
``` 
samples$Run

#[1] "SRR7457551_24_hours" "SRR7457552_24_hours" "SRR7457553_2_hours"  "SRR7457554_2_hours" 
#[5] "SRR7457555_2_hours"  "SRR7457556_2_hours"  "SRR7457557_control"  "SRR7457558_control" 
#[9] "SRR7457560_control"  "SRR7457562_24_hours" "SRR7457559_control"  "SRR7457561_24_hours"
``` 
### Set the path
``` 
file_path_stringtie <- file.path(dir,"stringtie_output", samples$Run)
``` 
### List files in the path
``` 
list.files(file_path_stringtie)

#[1] "e_data.ctab"                                 "e_data.ctab"                                
#[3] "e_data.ctab"                                 "e_data.ctab"                                
#[5] "e_data.ctab"                                 "e_data.ctab"                                
#[7] "e_data.ctab"                                 "e_data.ctab"                                
#[9] "e_data.ctab"                                 "e_data.ctab"                                
#[11] "e_data.ctab"                                 "e_data.ctab"                                
#[13] "e2t.ctab"                                    "e2t.ctab"       
``` 

### Assign the files from the stringtie_output
``` 
files_stringtie <- file.path(dir, "stringtie_output", samples$Run, "t_data.ctab")
files_stringtie

#[1] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/stringtie_output/SRR7457551_24_hours/t_data.ctab"
#[2] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/stringtie_output/SRR7457552_24_hours/t_data.ctab"
``` 
### Check if the files exist
``` 
file.exists(files_stringtie)

#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
``` 
### Assign sample names to files
``` 
names(files_stringtie) <- c('SRR7457551_24_hours','SRR7457552_24_hours','SRR7457553_2_hours','SRR7457554_2_hours','SRR7457555_2_hours','SRR7457556_2_hours','SRR7457557_control','SRR7457558_control','SRR7457560_control','SRR7457562_24_hours','SRR7457559_control','SRR7457561_24_hours')
``` 
### Read t_data.ctab as tmp
``` 
tmp <- read.table("/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/stringtie_output/SRR7457556_2_hours/t_data.ctab",header = TRUE)
``` 
### List head of tmp

head(tmp)

### Extract the two columns (t_name and gene_id) from data.frame contain transcripts ID and gene name information
``` 
#Assign tx2gene 
tx2gene <- tmp[, c("t_name", "gene_name")]

tx2gene

#                 t_name     gene_name
# 1    ENSMUST00000193812.1 4933401J01Rik
# 2    ENSMUST00000082908.1       Gm26206
# 3    ENSMUST00000162897.1          Xkr4
# 4    ENSMUST00000159265.1          Xkr4
# 5    ENSMUST00000070533.4          Xkr4
# 6    ENSMUST00000192857.1       Gm18956
# 7    ENSMUST00000195335.1       Gm37180
# 8    ENSMUST00000192336.1       Gm37363
``` 
### Import transcript-level abundance and counts for transcript

Different type of data generated by different software can be imported ("none", "salmon", "sailfish", "alevin", "kallisto", "rsem", "stringtie")

Here is the import type is stringtie
``` 
txi_stringtie <- tximport(files_stringtie, type = "stringtie",importer = NULL, tx2gene = tx2gene)

# reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9 10 11 12 
# summarizing abundance
# summarizing counts
# summarizing length

``` 
### Summary of the txi_stringtie
``` 
summary(txi_stringtie)

# Length Class  Mode     
# abundance           579852 -none- numeric  
# counts              579852 -none- numeric  
# length              579852 -none- numeric  
# countsFromAbundance      1 -none- character
``` 
### Take the counts data from the output from tximport
``` 
countData_stringtie <- txi_stringtie$counts
countData_stringtie

#SRR7457551_24_hours SRR7457552_24_hours SRR7457553_2_hours SRR7457554_2_hours
#0610006L08Rik           0.000000e+00        0.000000e+00       0.000000e+00       0.000000e+00
#0610007P14Rik           8.359200e+02        9.791066e+02       1.792400e+02       1.771733e+02
#0610009B22Rik           2.960533e+02        2.931133e+02       1.083333e+02       1.229733e+02
#0610009E02Rik           0.000000e+00        0.000000e+00       5.276726e+00       0.000000e+00
#0610009L18Rik           1.906667e+01        1.716000e+01       7.733332e+00       9.800000e+00
#0610009O20Rik           7.348667e+02        6.308400e+02       2.011602e+02       1.842728e+02
#0610010F05Rik           9.282868e+02        9.261133e+02       5.951334e+02       5.092534e+02
#0610010K14Rik           1.881463e+02        1.352396e+02       5.598808e+01       2.391762e+02
``` 
### Write the counts to a file
``` 
write.csv(countData_stringtie,file="countData_stringtie.csv")

``` 
### Analyzing RNA-Seq using DESeq2

[Go to Day2 Practical 5: Analyzing RNA-Seq using DESeq2](analyzing-RNA-seq-data-with-DESeq2.md)

## [Go to Day2 Practicals](rna-seq-wes-data-analysis-day2.md)

