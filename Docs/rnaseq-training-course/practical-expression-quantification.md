## StringTie Expression Quantification

> - StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. 

> - Output of StringTie can be processed by specialized software like Ballgown, Cuffdiff or other programs (DESeq2, edgeR, etc.).

> - In order to quantify transcript abundance from RNA-Seq data,  tools such as featureCounts, HTSeq are used.


``` 
#check the directory you are in by pwd command, it should be in the /Data_QC under your own scratch directory
[jl19@spectre13 Data_DC]$ pwd
/scratch/bbash/jl19/Data_QC
```
If it is already in your /Data_QC under your own scratch diretory e.g./scratch/bbash/jl19/Data_QC and then make "stringtie_output" directory, if not cd to the directory and then perform the command.
``` 
cd /scratch/bbash/jl19/Data_QC
``` 
``` 
[jl19@spectre13 Data_QC]$ mkdir stringtie_output
``` 

### Input files

#### Required files: sorted BAM and GFF files
> - Input files: sorted BAM files
> - General Feature Format (GFF)/General Transer Format(GTF):  Gene Annotation on the reference.

Both alignment files must be sorted by genomic location. The generic command line in this case becomes:

> - stringtie [-o <output.gtf>] --mix [other_options]  [short_read_alns.bam]  [long_read_alns.bam]

#### Option parameters for Stringtie

> - B option, it returns a Ballgown input table file, which contains coverage data for all transcripts.
> - Output file name (-o)
> - Reference annotation transcripts (-G)
> - Expression estimation mode (-e)

* When the -e option is used, the reference annotation file -G is a required input and StringTie will not attempt to assemble the input read alignments but instead it will only estimate the expression levels of the "reference" transcripts provided in the -G file.

* With this option, no "novel" transcript assemblies (isoforms) will be produced, and read alignments not overlapping any of the given reference transcripts will be ignored, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes for example.

> - gene_abund.tab (-A)	Gene abundances will be reported (tab delimited format) in the output file with the given name.


#### Download gencode.vM10.annotation.gff3.gz 

The file can be downloaded from  https://www.gencodegenes.org/mouse/release_M10.html

The file also can be copied from the server and unzip the file

```
cp -r  /data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/mm10/gencode.vM10.annotation.gff3.gz ./mm10

gunzip mm10/gencode.vM10.annotation.gff3.gz
```
---------

```
#stringtie_quantification.sh
##!/usr/bin/env bash

module load stringtie/2.1.1
#24 hours
stringtie trimmed/bam_files/SRR7457551_24_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457551_24_hours-trimmed.transcripts.gtf -A SRR7457551_24_hours-trimmed_gene_abund.tab;
stringtie trimmed/bam_files/SRR7457552_24_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457552_24_hours-trimmed.transcripts.gtf -A SRR7457551_24_hours-trimmed_gene_abund.tab;
stringtie trimmed/bam_files/SRR7457561_24_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457561_24_hours-trimmed.transcripts.gtf -A SRR7457551_24_hours-trimmed_gene_abund.tab;
stringtie SRR7457562_24_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457562_24_hours-trimmed.transcripts.gtf -A SRR7457562_24_hours-trimmed_gene_abund.tab;
#2 hours
stringtie SRR7457553_2_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457553_2_hours-trimmed.transcripts.gtf -A SRR7457553_2_hours-trimmed_gene_abund.tab;
stringtie SRR7457554_2_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457554_2_hours-trimmed.transcripts.gtf -A SRR7457554_2_hours-trimmed_gene_abund.tab;
stringtie SRR7457555_2_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457555_2_hours-trimmed.transcripts.gtf -A SRR7457555_2_hours-trimmed_gene_abund.tab;
stringtie SRR7457556_2_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457556_2_hours-trimmed.transcripts.gtf -A SRR7457556_2_hours-trimmed_gene_abund.tab;
#Control
stringtie SRR7457557_control-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457557_control-trimmed.transcripts.gtf -A SRR7457557_control-trimmed_gene_abund.tab;
stringtie SRR7457558_control-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457558_control-trimmed.transcripts.gtf -A SRR7457558_control-trimmed_gene_abund.tab;
stringtie SRR7457559_control-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457559_control-trimmed.transcripts.gtf -A SRR7457559_control-trimmed_gene_abund.tab;
stringtie SRR7457560_control-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457560_control-trimmed.transcripts.gtf -A SRR7457560_control-trimmed_gene_abund.tab;
```
### Output files

> - Generate gtf (containing the assembled transcripts that match the reference annotation) and t_data.ctab under the sample name folder

    stringtie_output/SRR7457551_24_hours/

    * SRR7457551_24_hours-trimmed.transcripts.gtf
    * e_data.ctab
    * e2t.ctab
    * i_data.ctab
    * i2t.ctab
    * t_data.ctab

> - stringtie_output/SRR7457552_24_hours-trimmed_gene_abund.tab

### Time Required to run each file

1.5min

### Prepare Sample list file as follow:
```
#Sample_list.txt

SRR7457551_24_hours	stringtie_output/SRR7457551_24_hours-trimmed.transcripts.gtf
SRR7457552_24_hours	stringtie_output/SRR7457552_24_hours-trimmed.transcripts.gtf
SRR7457561_24_hours	stringtie_output/SRR7457561_24_hours-trimmed.transcripts.gtf
SRR7457562_24_hours	stringtie_output/SRR7457562_24_hours-trimmed.transcripts.gtf

SRR7457553_2_hours	stringtie_output/SRR7457553_2_hours-trimmed.transcripts.gtf
SRR7457554_2_hours	stringtie_output/SRR7457554_2_hours-trimmed.transcripts.gtf
SRR7457555_2_hours	stringtie_output/SRR7457555_2_hours-trimmed.transcripts.gtf
SRR7457556_2_hours	stringtie_output/SRR7457556_2_hours-trimmed.transcripts.gtf

SRR7457557_control	stringtie_output/SRR7457557_control-trimmed.transcripts.gtf
SRR7457558_control	stringtie_output/SRR7457558_control-trimmed.transcripts.gtf
SRR7457559_control	stringtie_output/SRR7457559_control-trimmed.transcripts.gtf
SRR7457560_control	stringtie_output/SRR7457560_control-trimmed.transcripts.gtf
```

### Using Python script to process sample list to combine information from all samples
```
cp /data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/prepDE.py3  ./stringtie_output

cd stringtie_output
```
```
python prepDE.py3 Sample_list.txt
```

Two files are generated after running the script:
1. gene_count_matrix.csv

```
gene_id,SRR7457551_24_hours,SRR7457552_24_hours,SRR7457553_2_hours,SRR7457554_2_hours,SRR7457555_2_hours,SRR7457556_2_hours,SRR7457557_control,SRR7457558_control,SRR7457559_control,SRR7457560_control,SRR7457561_24_hours,SRR7457562_24_hours
ENSMUSG00000072660.5|Gm6288,4,0,8,2,0,0,0,0,0,0,0,0
ENSMUSG00000086111.1|Gm15326,0,0,2,0,0,0,0,0,0,0,0,0
ENSMUSG00000095516.1|Gm24135,1,0,0,0,0,0,1,0,1,1,1,0
ENSMUSG00000024386.8|Proc,0,0,0,0,0,0,0,0,0,0,0,0
ENSMUSG00000053310.11|Nrgn,6,6,5,0,7,0,0,0,0,9,7,12
ENSMUSG00000096150.1|Ighv1-85,0,0,0,0,0,0,0,0,0,0,0,0
ENSMUSG00000110019.1|RP23-453P1.2,0,0,0,0,1,1,1,1,0,0,0,0
ENSMUSG00000088229.1|Gm22604,0,0,0,0,0,0,0,0,0,0,0,0

.....
```
2. transcript_count_matrix.csv
```
transcript_id,SRR7457551_24_hours,SRR7457552_24_hours,SRR7457553_2_hours,SRR7457554_2_hours,SRR7457555_2_hours,SRR7457556_2_hours,SRR7457557_control,SRR7457558_control,SRR7457559_control,SRR7457560_control,SRR7457561_24_hours,SRR7457562_24_hours
ENSMUST00000119854.7,0,0,0,0,0,0,0,0,0,0,0,0
ENSMUST00000061309.4,0,0,0,0,0,0,0,0,0,0,0,0
ENSMUST00000070080.5,15,8,21,5,5,11,13,15,12,10,9,0
ENSMUST00000033908.13,11,11,18,10,22,32,11,5,13,16,24,12
ENSMUST00000174016.7,0,0,0,0,0,0,0,0,0,0,0,1
ENSMUST00000082868.1,0,0,0,0,0,0,0,0,0,0,0,0
ENSMUST00000187399.1,0,0,0,0,0,0,0,0,0,0,0,0
ENSMUST00000194793.1,2,4,5,0,7,5,8,6,4,8,3,3
ENSMUST00000072519.5,2329,1964,1494,1317,1891,1833,1208,1043,1471,1313,1760,1758
ENSMUST00000032272.12,23612,18864,19833,15290,19824,23624,15909,13396,19213,18874,19382,18569
...
```

### Extract Transcripts and Gene Information

The information can be found on public avalable databases: e.g. ensembl database.  Extract the two columns transcripts ID and gene name information.  
> Extract the two columns from the Mus_musculus.GRCm38.cdna.all.fa.gz which located at /scratch/bbash/jl19/Data_QC_May2023/mm10
> Extract all transcriptnames (1st) and genenames (4th) from  sequence names and write to a file.   

Here is short script to exact the above information from Mus_musculus.GRCm38.cdna.all.fa.gz:

```
#extract_transcriptname_genename.sh
gunzip -c Mus_musculus.GRCm38.cdna.all.fa.gz| \
grep '>' |
awk '{FS= " "}BEGIN{ print "TXNAME,GENEID"};{print substr($1,2) "," substr($4,6)};' > tx2gene.mm.GRCm38.cdna.csv 
```
```
cp -r /data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/mm10 ./
```
```
sh extract_transcriptname_genename.sh
```
Output file: tx2gene.mm.GRCm38.cdna.csv which contains TXNAME, GENEID

```
TXNAME,GENEID
ENSMUST00000177564.1,ENSMUSG00000096176.1
ENSMUST00000196221.1,ENSMUSG00000096749.2
ENSMUST00000179664.1,ENSMUSG00000096749.2
ENSMUST00000178537.1,ENSMUSG00000095668.1
ENSMUST00000178862.1,ENSMUSG00000094569.1
ENSMUST00000179520.1,ENSMUSG00000094028.1
ENSMUST00000179883.1,ENSMUSG00000094552.1
ENSMUST00000195858.1,ENSMUSG00000096420.2
ENSMUST00000179932.1,ENSMUSG00000096420.2
ENSMUST00000180001.1,ENSMUSG00000095656.1
```
### Tximport 

Import transcript-level estimates of StringTie and Kallisto output for DE analysis

[Go to Day2 Practical 4b: Tximport from StringTie and Kallisto Output](txtimport_StringTie_Kallisto.md)

### Analyzing RNA-Seq using DESeq2

[Go to Day2 Practical 5: Analyzing RNA-Seq using DESeq2](analyzing-RNA-seq-data-with-DESeq2.md)

## [Go to Day2 Practicals](rna-seq-wes-data-analysis-day2.md)
