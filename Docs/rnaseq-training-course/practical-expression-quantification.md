### StringTie Expression Quantification

> - StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. 

> - Output of StringTie can be processed by specialized software like Ballgown, Cuffdiff or other programs (DESeq2, edgeR, etc.).

> - In order to quantify transcript abundance from RNA-Seq data,  tools such as featureCounts, HTSeq are used.

``` 
module load stringtie/2.1.1
``` 
#### Input files

> - Required files: sorted BAM and GFF files
> - Input files: sorted BAM files
> - General Feature Format (GFF)/General Transer Format(GTF):  Gene Annotation on the reference.

Both alignment files must be sorted by genomic location. The generic command line in this case becomes:

> - stringtie [-o <output.gtf>] --mix [other_options]  [short_read_alns.bam]  [long_read_alns.bam]

The regular options include a reference annotation (-G) and an output file name (-o), so a more realistic command line example would look like this:

> - Download gencode.vM10.annotation.gff.gz from  https://www.gencodegenes.org/mouse/release_M10.html

```
gunzip gencode.vM10.annotation.gff.gz
```

The file also can be copied from 

```
cp -r  /data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/mm10/gencode.vM10.annotation.gff.gz ./

```
> - Reference annotation transcripts (-G)

> - Expression estimation mode (-e)

When the -e option is used, the reference annotation file -G is a required input and StringTie will not attempt to assemble the input read alignments but instead it will only estimate the expression levels of the "reference" transcripts provided in the -G file.

With this option, no "novel" transcript assemblies (isoforms) will be produced, and read alignments not overlapping any of the given reference transcripts will be ignored, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes for example.

---------
>  Generate gtf and t_data.ctab under the sample name folder

>  stringtie_output/SRR7457551/

> - e_data.ctab
> - e2t.ctab
> - i_data.ctab
> - i2t.ctab
> - t_data.ctab

```
stringtie SRR7457551.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457551.transcripts.gtf -A SRR7457551_gene_abund.tab;
stringtie SRR7457552.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457552.transcripts.gtf -A SRR7457552_gene_abund.tab;
stringtie SRR7457555.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457555.transcripts.gtf -A SRR7457555_gene_abund.tab;
stringtie SRR7457556.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457556.transcripts.gtf -A SRR7457556_gene_abund.tab;
stringtie SRR7457559.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457559.transcripts.gtf -A SRR7457559_gene_abund.tab;
stringtie SRR7457560.sorted.bam -B -G  mm10/gencode.vM10.annotation.gff3 -e -o stringtie_output/SRR7457560.transcripts.gtf -A SRR7457560_gene_abund.tab;

```
#### Output files

> - GTF file: containing the assembled transcripts that match the reference annotation

### Prepare Sample list file as follow:
>  Sample_list.txt

> - SRR7457551	stringtie_output/SRR7457551.transcripts.gtf
> - SRR7457552	stringtie_output/SRR7457552.transcripts.gtf
> - SRR7457555	stringtie_output/SRR7457555.transcripts.gtf
> - SRR7457556	stringtie_output/SRR7457556.transcripts.gtf
> - SRR7457559	stringtie_output/SRR7457559.transcripts.gtf
> - SRR7457560	stringtie_output/SRR7457560.transcripts.gtf

#### Using Python script to process sample list
```
python prepDE.py3 stringtie_output/sample_list.txt
```
https://nasqar.abudhabi.nyu.edu/deseq2shiny/


#### Quantification

Extract the two columns transcripts ID and gene name information, The information can be found on public avalable databases: e.g. ensembl database

> Extract the two columns from the Mus_musculus.GRCm39.cdna.all.fa.gz

> Extract all transcriptnames (1st) and genenames (4th) from  sequence names and write to a file.   
[extract_script](extract_script.sh)


#### Tximport 

Import transcript-level estimates of StringTie and Kallisto output for DE analysis
[Tximport from StringTie and Kallisto Output](txtimport_StringTie_Kallisto.md)


## [Go to Day2 Practicals](rna-seq-wes-data-analysis-day2#quantification)
