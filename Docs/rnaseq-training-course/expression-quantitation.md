#Expression Quantition

StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts. 

StringTie's output can be processed by specialized software like Ballgown, Cuffdiff or other programs (DESeq2, edgeR, etc.).
St
In order to quantify transcript abundance from RNA-Seq data, we need to use tools featureCpunts, HTSeq  

``` 
module load stringtie/2.1.1
``` 
###Input files

sorted BAM files can be used as input 



Required files: BAM and GFF file

Both alignment files must be sorted by genomic location. The generic command line in this case becomes:

stringtie [-o <output.gtf>] --mix [other_options] <short_read_alns.bam> <long_read_alns.bam>
The regular options include a reference annotation (-G) and an output file name (-o), so a more realistic command line example would look like this:

Download gencode.vM10.annotation.gff.gz from  https://www.gencodegenes.org/mouse/release_M10.html
```
gunzip gencode.vM10.annotation.gff.gz
```

```
stringtie SRR7457551.sorted.bam -G mm10/gencode.vM10.annotation.gff3  -o SRR7457551.transcripts.gtf;
stringtie SRR7457552.sorted.bam -G mm10/gencode.vM10.annotation.gff3  -o SRR7457552.transcripts.gtf;
stringtie SRR7457555.sorted.bam -G mm10/gencode.vM10.annotation.gff3  -o SRR7457555.transcripts.gtf;
stringtie SRR7457556.sorted.bam -G mm10/gencode.vM10.annotation.gff3  -o SRR7457556.transcripts.gtf;
stringtie SRR7457559.sorted.bam -G mm10/gencode.vM10.annotation.gff3  -o SRR7457559.transcripts.gtf;
stringtie SRR7457560.sorted.bam -G mm10/gencode.vM10.annotation.gff3  -o SRR7457560.transcripts.gtf;
```

 