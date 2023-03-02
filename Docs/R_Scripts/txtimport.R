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

BiocManager::install("tximport")
BiocManager::install("haven")
#BiocManager::install("problems")
library("tximport")

#set directory

dir <- "~/BASH_NGS_Course/docs/RNA-seq-example"

#list files in the directory

list.files(dir)

#set the path

file_path <- file.path(dir, "stringtie_output")

#list files in the path
#list.files(path = "~/BASH_NGS_Course/docs/RNA-seq-example/stringtie_output")
list.files(file_path)


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

#rownames(samples) <- samples$Run

#Set the factor as 24h, 2h and Naive replicates 2
samples$condition <- factor(rep(c("24h","2h","Naive"),each =2))

# [1] 24h   24h   2h    2h    Naive Naive
# Levels: 24h 2h Naive

#Read t_data.ctab as tmp

tmp <- read.table("~/BASH_NGS_Course/docs/RNA-seq-example/stringtie_output/SRR7457555/t_data.ctab",header = TRUE)

tmp
# t_id  chr strand   start     end                t_name num_exons length               gene_id     gene_name       cov      FPKM
# 1     1 chr1      + 3073253 3074322  ENSMUST00000193812.1         1   1070  ENSMUSG00000102693.1 4933401J01Rik  0.000000  0.000000
# 2     2 chr1      + 3102016 3102125  ENSMUST00000082908.1         1    110  ENSMUSG00000064842.1       Gm26206  0.000000  0.000000
# 3     3 chr1      - 3205901 3216344  ENSMUST00000162897.1         2   4153  ENSMUSG00000051951.5          Xkr4  0.000000  0.000000
# 4     4 chr1      - 3206523 3215632  ENSMUST00000159265.1         2   2989  ENSMUSG00000051951.5          Xkr4  0.000000  0.000000
# 5     5 chr1      - 3214482 3671498  ENSMUST00000070533.4         3   3634  ENSMUSG00000051951.5          Xkr4  0.000000  0.000000
# 6     6 chr1      + 3252757 3253236  ENSMUST00000192857.1         1    480  ENSMUSG00000102851.1       Gm18956  0.000000  0.000000
# 7     7 chr1      - 3365731 3368549  ENSMUST00000195335.1         1   2819  ENSMUSG00000103377.1       Gm37180  0.000000  0.000000

#read all the t_data.ctab files under the diretory: ~/BASH_NGS_Course/docs/RNA-seq-example/stringtie_output/SRR7457551/t_data.ctab
#There are t_data.ctab in each direcotry with different sample names:SRR7457551','SRR7457552','SRR7457555','SRR7457556','SRR7457559','SRR7457560'
files <- file.path(dir, "stringtie_output", samples$Run, "t_data.ctab")

#Check if the files exist
file.exists(files)

#[1] TRUE TRUE TRUE TRUE TRUE TRUE

#Assign sample names to files
names(files) <- c('SRR7457551','SRR7457552','SRR7457555','SRR7457556','SRR7457559','SRR7457560')

#list head of tmp
head(tmp)

#Extract the two columns  data.frame contain transcripts ID and gene name information: t_name, gene_name

tx2gene <- tmp[, c("t_name", "gene_name")]

tx2gene

# t_name     gene_name
# 1    ENSMUST00000193812.1 4933401J01Rik
# 2    ENSMUST00000082908.1       Gm26206
# 3    ENSMUST00000162897.1          Xkr4
# 4    ENSMUST00000159265.1          Xkr4
# 5    ENSMUST00000070533.4          Xkr4
# 6    ENSMUST00000192857.1       Gm18956
# 7    ENSMUST00000195335.1       Gm37180
# 8    ENSMUST00000192336.1       Gm37363

#Import transcript-level abundance and counts for transcript
#Different type of data generated by different software can be imported ("none", "salmon", "sailfish", "alevin", "kallisto", "rsem", "stringtie")

#Here is the import type is stringtie
txi <- tximport(files, type = "stringtie",importer = NULL, tx2gene = tx2gene)

# reading in files with read_tsv
# 1 2 3 4 5 6 
# summarizing abundance
# summarizing counts
# summarizing length


#Summary of the txi
summary(txi)

# Length Class  Mode     
# abundance           289926 -none- numeric  
# counts              289926 -none- numeric  
# length              289926 -none- numeric  
# countsFromAbundance      1 -none- character

#Take the counts data from the output from tximport
countData <- txi$counts

#Write the counts to a file
write.csv(countData,file="countData_txiimport.csv")

