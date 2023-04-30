### StringTie Expression Quantification

> - StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. 

> - Output of StringTie can be processed by specialized software like Ballgown, Cuffdiff or other programs (DESeq2, edgeR, etc.).

> - In order to quantify transcript abundance from RNA-Seq data,  tools such as featureCounts, HTSeq are used.


``` 

#check the directory you are in by pwd command
[jl19@spectre13 Data_DC]$ pwd
/scratch/bbash/jl19/Data_QC
```
#if it is already in /scratch/bbash/jl19/Data_QC and make "stringtie_output" directory, if not 
``` 
cd /scratch/bbash/jl19/Data_QC
``` 
``` 
[jl19@spectre13 Data_QC]$ mkdir stringtie_output

``` 

#### Input files

> - Required files: sorted BAM and GFF files
> - Input files: sorted BAM files
> - General Feature Format (GFF)/General Transer Format(GTF):  Gene Annotation on the reference.

Both alignment files must be sorted by genomic location. The generic command line in this case becomes:

> - stringtie [-o <output.gtf>] --mix [other_options]  [short_read_alns.bam]  [long_read_alns.bam]

The regular options include a reference annotation (-G) and an output file name (-o), so a more realistic command line example would look like this:

> - Download gencode.vM10.annotation.gff3.gz from  https://www.gencodegenes.org/mouse/release_M10.html


The file also can be copied from the server and unzip the file

```
cp -r  /data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/mm10/gencode.vM10.annotation.gff3.gz ./mm10

gunzip mm10/gencode.vM10.annotation.gff3.gz

```
> - Reference annotation transcripts (-G)

> - Expression estimation mode (-e)

When the -e option is used, the reference annotation file -G is a required input and StringTie will not attempt to assemble the input read alignments but instead it will only estimate the expression levels of the "reference" transcripts provided in the -G file.

With this option, no "novel" transcript assemblies (isoforms) will be produced, and read alignments not overlapping any of the given reference transcripts will be ignored, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes for example.

> - gene_abund.tab (-A)	Gene abundances will be reported (tab delimited format) in the output file with the given name.

---------

```
#stringtie_quantification.sh
##!/usr/bin/env bash

module load stringtie/2.1.1

stringtie SRR7457551_24_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457551_24_hours-trimmed.transcripts.gtf -A SRR7457551_24_hours-trimmed_gene_abund.tab;
stringtie SRR7457552_24_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457552_24_hours-trimmed.transcripts.gtf -A SRR7457551_24_hours-trimmed_gene_abund.tab;
stringtie SRR7457561_24_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457561_24_hours-trimmed.transcripts.gtf -A SRR7457551_24_hours-trimmed_gene_abund.tab;
stringtie SRR7457562_24_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457562_24_hours-trimmed.transcripts.gtf -A SRR7457562_24_hours-trimmed_gene_abund.tab;


stringtie SRR7457553_2_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457553_2_hours-trimmed.transcripts.gtf -A SRR7457553_2_hours-trimmed_gene_abund.tab;
stringtie SRR7457554_2_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457554_2_hours-trimmed.transcripts.gtf -A SRR7457554_2_hours-trimmed_gene_abund.tab;
stringtie SRR7457555_2_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457555_2_hours-trimmed.transcripts.gtf -A SRR7457555_2_hours-trimmed_gene_abund.tab;
stringtie SRR7457556_2_hours-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457556_2_hours-trimmed.transcripts.gtf -A SRR7457556_2_hours-trimmed_gene_abund.tab;

stringtie SRR7457557_control-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457557_control-trimmed.transcripts.gtf -A SRR7457557_control-trimmed_gene_abund.tab;
stringtie SRR7457558_control-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457558_control-trimmed.transcripts.gtf -A SRR7457558_control-trimmed_gene_abund.tab;
stringtie SRR7457559_control-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457559_control-trimmed.transcripts.gtf -A SRR7457559_control-trimmed_gene_abund.tab;
stringtie SRR7457560_control-trimmed.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3.gz -e -o stringtie_output/SRR7457560_control-trimmed.transcripts.gtf -A SRR7457560_control-trimmed_gene_abund.tab;


```
#### Output files

1. Generate gtf (containing the assembled transcripts that match the reference annotation) and t_data.ctab under the sample name folder

    >  stringtie_output/SRR7457551_24_hours/

    > - SRR7457551_24_hours-trimmed.transcripts.gtf
    > - e_data.ctab
    > - e2t.ctab
    > - i_data.ctab
    > - i2t.ctab
    > - t_data.ctab

2. stringtie_output/SRR7457552_24_hours-trimmed_gene_abund.tab

### Time Required to run each file

1.5min

### Prepare Sample list file as follow:
>  Sample_list.txt

> - SRR7457551_24_hours	stringtie_output/SRR7457551_24_hours-trimmed.transcripts.gtf
> - SRR7457552_24_hours	stringtie_output/SRR7457552_24_hours-trimmed.transcripts.gtf
> - SRR7457561_24_hours	stringtie_output/SRR7457561_24_hours-trimmed.transcripts.gtf
> - SRR7457562_24_hours	stringtie_output/SRR7457562_24_hours-trimmed.transcripts.gtf

> - SRR7457553_2_hours	stringtie_output/SRR7457553_2_hours-trimmed.transcripts.gtf
> - SRR7457554_2_hours	stringtie_output/SRR7457554_2_hours-trimmed.transcripts.gtf
> - SRR7457555_2_hours	stringtie_output/SRR7457555_2_hours-trimmed.transcripts.gtf
> - SRR7457556_2_hours	stringtie_output/SRR7457556_2_hours-trimmed.transcripts.gtf

> - SRR7457557_control	stringtie_output/SRR7457557_control-trimmed.transcripts.gtf
> - SRR7457558_control	stringtie_output/SRR7457558_control-trimmed.transcripts.gtf

> - SRR7457559_control	stringtie_output/SRR7457559_control-trimmed.transcripts.gtf
> - SRR7457560_control	stringtie_output/SRR7457560_control-trimmed.transcripts.gtf


#### Using Python script to process sample list
```
cp /data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/prepDE.pye  ./stringtie_output

cd stringtie_output
```
```
python prepDE.py3 sample_list.txt
```

Two files are generated:
1. gene_count_matrix.csv
2. transcript_count_matrix.csv

https://nasqar.abudhabi.nyu.edu/deseq2shiny/


#### Quantification

Extract the two columns transcripts ID and gene name information.  The information can be found on public avalable databases: e.g. ensembl database

> Extract the two columns from the Mus_musculus.GRCm39.cdna.all.fa.gz which located at /scratch/bbash/jl19/Data_QC_May2023/mm10
> Extract all transcriptnames (1st) and genenames (4th) from  sequence names and write to a file.   
[extract_script](extract_script.sh)


#### Tximport 

Import transcript-level estimates of StringTie and Kallisto output for DE analysis
[Tximport from StringTie and Kallisto Output](txtimport_StringTie_Kallisto.md)

[Go to Day2 Practical 5](analyzing-RNA-seq-data-with-DESeq2.md)
## [Go to Day2 Practicals](rna-seq-wes-data-analysis-day2.md)
