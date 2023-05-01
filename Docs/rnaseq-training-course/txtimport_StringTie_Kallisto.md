## Import Transcript-level Estimates

> - Imports transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages. 
> - Average transcript length, weighted by sample-specific transcript abundance estimates, is provided as a matrix which can be used as an offset for different expression of gene-level counts.


``` 
#################################
##Import from Kallisto Output##
#################################

#Install required R packages
if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install("ballgown")
BiocManager::install("rhdf5")
BiocManager::install("tximport")
BiocManager::install("haven")

#################################
##Import from Kallisto Output##
#################################


library(dplyr) # data wrangling
library(ggplot2) # plotting
library(DESeq2) # rna-seq
library(edgeR) # rna-seq
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5) # read/convert kalisto output files.  



#################################
##Import from kallisto Output##
#################################

#set working directory
dir <- "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583"

setwd("/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/DE_Analysis/")


#Read sample information
samples <- read.csv("sample_metadata_12samples_kallisto.csv", header = TRUE)
#sample run information
samples$Run

#[1] "SRR7457551_24_hours-trimmed" "SRR7457552_24_hours-trimmed" "SRR7457553_2_hours-trimmed" 
#[4] "SRR7457554_2_hours-trimmed"  "SRR7457555_2_hours-trimmed"  "SRR7457556_2_hours-trimmed" 
#[7] "SRR7457557_control-trimmed"  "SRR7457558_control-trimmed"  "SRR7457560_control-trimmed" 
#[10] "SRR7457562_24_hours-trimmed" "SRR7457559_control-trimmed"  "SRR7457561_24_hours-trimmed"
#set the path

file_path_kallisto <- file.path(dir,"kallisto_output", samples$Run)

#list files in the path
list.files(file_path_kallisto)

#[1] "abundance.h5"  "abundance.h5"  "abundance.h5"  "abundance.h5"  "abundance.h5"  "abundance.h5" 
#[7] "abundance.h5"  "abundance.h5"  "abundance.h5"  "abundance.h5"  "abundance.h5"  "abundance.tsv"
#[13] "abundance.tsv" "abundance.tsv" "abundance.tsv" "abundance.tsv" "abundance.tsv" "abundance.tsv"
#[19] "abundance.tsv" "abundance.tsv" "abundance.tsv" "abundance.tsv" "run_info.json" "run_info.json"
#[25] "run_info.json" "run_info.json" "run_info.json" "run_info.json" "run_info.json" "run_info.json"
#[31] "run_info.json" "run_info.json" "run_info.json"

files_kallisto <- file.path(dir, "kallisto_output", samples$Run, "abundance.tsv")
files_kallisto
# [1] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583//kallisto_output/SRR7457551_24_hours-trimmed/abundance.tsv"
# [2] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583//kallisto_output/SRR7457552_24_hours-trimmed/abundance.tsv"
#....
# [11] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457559_control-trimmed/abundance.tsv" 
# [12] "/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/SRR7457561_24_hours-trimmed/abundance.tsv"


#Check if the files exist
file.exists(files_kallisto)

#[1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

#Assign sample names to files
names(files_kallisto) <- c('SRR7457551_24_hours','SRR7457552_24_hours','SRR7457553_2_hours','SRR7457554_2_hours','SRR7457555_2_hours','SRR7457556_2_hours','SRR7457557_control','SRR7457558_control','SRR7457560_control','SRR7457562_24_hours','SRR7457559_control','SRR7457561_24_hours')



#List the sample information
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
     
#List the run names

samples$Run

#[1] "SRR7457551_24_hours-trimmed" "SRR7457552_24_hours-trimmed" "SRR7457553_2_hours-trimmed" 
#[4] "SRR7457554_2_hours-trimmed"  "SRR7457555_2_hours-trimmed"  "SRR7457556_2_hours-trimmed" 
#[7] "SRR7457557_control-trimmed"  "SRR7457558_control-trimmed"  "SRR7457560_control-trimmed" 
#[10] "SRR7457562_24_hours-trimmed" "SRR7457559_control-trimmed"  "SRR7457561_24_hours-trimmed"




#Transcripts need to be associated with gene IDs for gene-level summarization. If that information is present in the files, we can skip this step. For Salmon, Sailfish, and kallisto the files only provide the transcript ID. We first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. The column names do not matter but this column order must be used. The transcript ID must be the same one used in the abundance files.

#Creating this tx2gene data.frame can be accomplished from a TxDb object and the select function from the AnnotationDbi package. The following code could be used to construct such a table:
  

#####################################################################################################################################
#For kallisto the files only provide the transcript ID. We first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. The column names do not matter but this column order must be used. The transcript ID must be the same one used in the abundance files.

#Creating this tx2gene data.frame can be accomplished from a TxDb object and the select function from the AnnotationDbi package. The following code could be used to construct such a table:
#####################################################################################################################################################

#Load TxDb annotation package

BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", force = TRUE)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
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


columns(txdb)
# [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSPHASE"   "CDSSTART"   "CDSSTRAND"  "EXONCHROM"  "EXONEND"    "EXONID"     "EXONNAME"  
# [12] "EXONRANK"   "EXONSTART"  "EXONSTRAND" "GENEID"     "TXCHROM"    "TXEND"      "TXID"       "TXNAME"     "TXSTART"    "TXSTRAND"   "TXTYPE"    

keytypes(txdb) 

## [1] "CDSID"    "CDSNAME"  "EXONID"   "EXONNAME" "GENEID"   "TXID"     "TXNAME"  


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
countData <- txi_kallisto$counts

countData
# SRR7457551_24_hours SRR7457552_24_hours SRR7457553_2_hours SRR7457554_2_hours
                  9.899398e+05        9.752320e+05       8.422768e+05       6.823460e+05
# 100009600        1.903941e+01        2.174964e+01       1.623218e+01       2.352240e+00
# 100009609        1.138890e+01        4.314980e+00       7.382950e+00       7.280110e+00
# 100009614        0.000000e+00        0.000000e+00       0.000000e+00       0.000000e+00
# 100012           0.000000e+00        0.000000e+00       0.000000e+00       0.000000e+00
# 100017           1.495003e+03        1.263000e+03       1.720005e+03       1.088000e+03
# 100019           2.076546e+03        2.019474e+03       1.960450e+03       1.199829e+03

#Write the counts to a file
write.csv(countData,file="countData_kallisto.csv")


``` 
