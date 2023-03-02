#########################################
####Import transcript-level estimates####
#########################################


#Imports transcript-level abundance, estimated counts and transcript lengths, and 
#summarizes into matrices for use with downstream gene-level analysis packages. 
#Average transcript length, weighted by sample-specific transcript abundance estimates, 
#is provided as a matrix which can be used as an offset for different expression of gene-level counts.

if (!requireNamespace("BiocManager", quietly=TaRUE))
install.packages("BiocManager")
BiocManager::install("ballgown")
BiocManager::install("rhdf5")
BiocManager::install("tximport")
BiocManager::install("haven")

#################################
##Improt from Kallisto Output##
#################################


library(dplyr) # data wrangling
library(ggplot2) # plotting
library(DESeq2) # rna-seq
library(edgeR) # rna-seq
library(tximport) # importing kalisto transcript counts to geneLevels
library(readr) # Fast readr of files.
library(rhdf5) # read/convert kalisto output files.  



#################################
##Improt from lallisto Output##
#################################

#set directory

dir <- "~/BASH_NGS_Course/docs/RNA-seq-example"

#set the path

file_path_kallisto <- file.path(dir, "kallisto_output", samples$Run)
#list files in the path

list.files(file_path_kallisto)

files_kallisto <- file.path(dir, "kallisto_output", samples$Run, "abundance.tsv")
files_kallisto
# [1] "~/BASH_NGS_Course/docs/RNA-seq-example/kallisto_output/SRR7457551/abundance.tsv" "~/BASH_NGS_Course/docs/RNA-seq-example/kallisto_output/SRR7457552/abundance.tsv"
# [3] "~/BASH_NGS_Course/docs/RNA-seq-example/kallisto_output/SRR7457555/abundance.tsv" "~/BASH_NGS_Course/docs/RNA-seq-example/kallisto_output/SRR7457556/abundance.tsv"
# [5] "~/BASH_NGS_Course/docs/RNA-seq-example/kallisto_output/SRR7457559/abundance.tsv" "~/BASH_NGS_Course/docs/RNA-seq-example/kallisto_output/SRR7457560/abundance.tsv"

#Check if the files exist
file.exists(files_kallisto)

#[1] TRUE TRUE TRUE TRUE TRUE TRUE

#Assign sample names to files
names(files_kallisto) <- c('SRR7457551','SRR7457552','SRR7457555','SRR7457556','SRR7457559','SRR7457560')


#Read sample information
samples <- read.csv("~/BASH_NGS_Course/docs/RNA-seq-example/sample_information.csv", header = TRUE)

#List the sample information
samples

# Run      Age Assay.Type AvgSpotLen      Bases  BioProject    BioSample           BioSampleModel     Bytes           Cell_type             Center.Name Consent DATASTORE.filetype
# 1 SRR7457551 14 weeks    RNA-Seq         75 2244859888 PRJNA450151 SAMN08949698 Model organism or animal 883354260 Alveolar Macrophage NORTHWESTERN UNIVERSITY  public          sra,fastq
# 2 SRR7457552 14 weeks    RNA-Seq         75 2018541831 PRJNA450151 SAMN08949699 Model organism or animal 793650717 Alveolar Macrophage NORTHWESTERN UNIVERSITY  public          sra,fastq
# 3 SRR7457555 14 weeks    RNA-Seq         75 2011980926 PRJNA450151 SAMN08949694 Model organism or animal 788409431 Alveolar Macrophage NORTHWESTERN UNIVERSITY  public          fastq,sra
# 4 SRR7457556 14 weeks    RNA-Seq         75 2029890687 PRJNA450151 SAMN08949695 Model organism or animal 796029883 Alveolar Macrophage NORTHWESTERN UNIVERSITY  public          sra,fastq
# 5 SRR7457559 14 weeks    RNA-Seq         75 1703298074 PRJNA450151 SAMN08949690 Model organism or animal 666357446 Alveolar Macrophage NORTHWESTERN UNIVERSITY  public          sra,fastq
# 6 SRR7457560 14 weeks    RNA-Seq         75 1511125335 PRJNA450151 SAMN08949691 Model organism or animal 591158034 Alveolar Macrophage NORTHWESTERN UNIVERSITY  public          sra,fastq
# DATASTORE.provider   DATASTORE.region Experiment  Instrument                                   Isolate   Library.Name LibraryLayout LibrarySelection  LibrarySource     Organism Platform
# 1              s3,gs s3.us-east-1,gs.US SRX4328058 NextSeq 500 24h Reperfusion AM Biological Replicate 1 R05_AM_Allo24h        SINGLE              PCR TRANSCRIPTOMIC Mus musculus ILLUMINA
# 2              gs,s3 s3.us-east-1,gs.US SRX4328057 NextSeq 500 24h Reperfusion AM Biological Replicate 2 R06_AM_Allo24h        SINGLE              PCR TRANSCRIPTOMIC Mus musculus ILLUMINA
# 3              gs,s3 s3.us-east-1,gs.US SRX4328054 NextSeq 500  2h Reperfusion AM Biological Replicate 1  R01_AM_Allo2h        SINGLE              PCR TRANSCRIPTOMIC Mus musculus ILLUMINA
# 4              s3,gs gs.US,s3.us-east-1 SRX4328053 NextSeq 500  2h Reperfusion AM Biological Replicate 2  R02_AM_Allo2h        SINGLE              PCR TRANSCRIPTOMIC Mus musculus ILLUMINA
# 5              gs,s3 s3.us-east-1,gs.US SRX4328050 NextSeq 500           Naive AM Biological Replicate 1   N01_AM_Naive        SINGLE              PCR TRANSCRIPTOMIC Mus musculus ILLUMINA
# 6              s3,gs gs.US,s3.us-east-1 SRX4328049 NextSeq 500           Naive AM Biological Replicate 2   N02_AM_Naive        SINGLE              PCR TRANSCRIPTOMIC Mus musculus ILLUMINA
# ReleaseDate    Sample.Name  sex SRA.Study        strain tissue
# 1 2018-06-29T00:00:00Z R05_AM_Allo24h male SRP151689 Cx3cr1gfp/+B6   Lung
# 2 2018-06-29T00:00:00Z R06_AM_Allo24h male SRP151689 Cx3cr1gfp/+B6   Lung
# 3 2018-06-29T00:00:00Z  R01_AM_Allo2h male SRP151689 Cx3cr1gfp/+B6   Lung
# 4 2018-06-29T00:00:00Z  R02_AM_Allo2h male SRP151689 Cx3cr1gfp/+B6   Lung
# 5 2018-06-29T00:00:00Z   N01_AM_Naive male SRP151689 Cx3cr1gfp/+B6   Lung
# 6 2018-06-29T00:00:00Z   N02_AM_Naive male SRP151689 Cx3cr1gfp/+B6   Lung
# 
#List the run names

samples$Run

#[1] "SRR7457551" "SRR7457552" "SRR7457555" "SRR7457556" "SRR7457559" "SRR7457560"

#Set the factor as 24h, 2h and Naive replicates 2
samples$condition <- factor(rep(c("24h","2h","Naive"),each =2))

# [1] 24h   24h   2h    2h    Naive Naive
# Levels: 24h 2h Naive

#Read t_data.ctab as tmp

tmp <- read.table("~/BASH_NGS_Course/docs/RNA-seq-example/kallisto_output/SRR7457555/abundance.tsv",header = TRUE)

tmp
# target_id length eff_length  est_counts         tpm
# 1    ENSMUST00000178537.2     12    4.16231 0.00000e+00 0.00000e+00
# 2    ENSMUST00000178862.2     14    4.39563 0.00000e+00 0.00000e+00
# 3    ENSMUST00000196221.2      9    3.70903 0.00000e+00 0.00000e+00
# 4    ENSMUST00000179664.2     11    4.02702 0.00000e+00 0.00000e+00
#read all the abundance.tsv files under the diretory: ~/BASH_NGS_Course/docs/RNA-seq-example/kallisto_output/SRR7457551/abundance.tsv
#There are t_data.ctab in each direcotry with different sample names:SRR7457551','SRR7457552','SRR7457555','SRR7457556','SRR7457559','SRR7457560'
#[1] TRUE TRUE TRUE TRUE TRUE TRUE


#list head of tmp
head(tmp)


#Transcripts need to be associated with gene IDs for gene-level summarization. If that information is present in the files, we can skip this step. For Salmon, Sailfish, and kallisto the files only provide the transcript ID. We first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. The column names do not matter but this column order must be used. The transcript ID must be the same one used in the abundance files.

#Creating this tx2gene data.frame can be accomplished from a TxDb object and the select function from the AnnotationDbi package. The following code could be used to construct such a table:
  

#Extract the two columns  data.frame contain transcripts ID and gene name information: t_name, gene_name

# Extract all transcriptnames (1st) and genenames (4th) from  sequence names and write to a file.   


###########################################################################################
# gzcat ~/BASH_NGS_Course/docs/RNA-seq-example/Mus_musculus.GRCm39.cdna.all.fa.gz| \
# grep '>' |
#   awk '{FS= " "}BEGIN{ print "t_name,gene_name"};{print substr($1,2) "," substr($7,13)};' \
# > ~/BASH_NGS_Course/docs/RNA-seq-example/tx2gene.mm.GRCm39.cdna.csv 
# 
###########################################################################################


#####################################################################################################################################
#For Salmon, Sailfish, and kallisto the files only provide the transcript ID. We first make a data.frame called tx2gene with two columns: 1) transcript ID and 2) gene ID. The column names do not matter but this column order must be used. The transcript ID must be the same one used in the abundance files.

#Creating this tx2gene data.frame can be accomplished from a TxDb object and the select function from the AnnotationDbi package. The following code could be used to construct such a table:
######################################################################################################################################################################################

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

#samples <- read.csv("~/BASH_NGS_Course/docs/RNA-seq-example/sample_information.csv",header = TRUE)

tx2gene_db

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
txi_kallisto <- tximport(files, type = "kallisto", tx2gene = tx2gene_db, ignoreAfterBar = TRUE,dropInfReps = TRUE)

# 1 2 3 4 5 6 
# transcripts missing from tx2gene: 117936
# summarizing abundance
# summarizing counts
# summarizing length

#Summary of the txi
summary(txi_kallisto)

# Length Class  Mode     
# abundance           78 -none- numeric  
# counts              78 -none- numeric  
# length              78 -none- numeric  
# countsFromAbundance      1 -none- character

#Take the counts data from the output from tximport
countData <- txi_kallisto$counts

#Write the counts to a file
write.csv(countData,file="countData_kallisto.csv")


#tx2gene information can also be imported directly from file:gencode.vM18.annotation.tx2gene
#file:gencode.vM18.annotation.tx2gene.csv can be generated by a R script or download diretly from
#


