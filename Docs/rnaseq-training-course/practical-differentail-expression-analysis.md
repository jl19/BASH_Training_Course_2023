##Expression levels, Read Count and Normalisation
  
###Raw read counts are not informative due to different sources of bias: gene length, GC-content and sequencing depth;
  
Most common normatlization metrics;  
  
| Abbreviation|Meaning         | Explanation  |                      
| ----------- | ---------------------|--------------|
| `RPKMs`     | reads per kilobase of transript per million reads mapped  | a normalized unit of transcripts/genes expression, it scales by transcriopt length |
| `pm`| per million scale factor|[total number of reads]/1,000,000 |
| `Gene A RPM`       | =[Gene A read counts]/pm |This normalizes for sequencing depth, giving you reads per million (RPM)               |
| `Gene A RPKM`     |[Gene A RPM]/[length of Gene A]      | Divide the RPM values by the length of the gene, in kilobases.This gives you RPKM.|
|`FPKMs` |fragments per kilobase of transcript per million reads mapped |FPKM is very similar to RPKM. RPKM was made for single-end RNA-seq, where every read corresponded to a single fragment that was sequenced. 
FPKM was made for paired-end RNA-seq. With paired-end RNA-seq, two reads can correspond to a single fragment, or, if one read in the pair did not map,one read can correspond to a single fragment. The only difference between RPKM and FPKM is that FPKM takes into account that two reads can map to one fragment|
|`TPM`|(transcripts per million)|TPM is very similar to RPKM and FPKM. The only difference is the order of operations. |
|`Gene A RPK` |= [Gene A raw read counts]/[length of gene A] |Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).|
pm = RPK/1,000,000. Count up all the RPK values in a sample and divide this number by 1,000,000. This is  “per million” scaling factor.
|`Gene A TPM`|= Gene A RPK/pm. |Divide the RPK values by the “per million” scaling factor. This gives you TPM.||
|When you use TPM, the sum of all TPMs in each sample are the same. This makes it easier to compare the proportion of reads that mapped to a gene in each sample. |
|In contrast, with RPKM and FPKM, the sum of the normalized reads in each sample may be different, and this makes it harder to compare samples directly.|
|DESeq2 using RLE(Relative Log Expression) and edgeR using TMM(Trimmed Mean of M-values) normalization methods for analysing RNA-seq data, they were shown to produce quantifications robust to the presence of different library sizes and widely different libary compositions.|


##StringTie

use to normalize the data

##Differential Gene Expression Analysis

[DESeq2](analyzing-RNA-seq-data-with-DESeq2.md)
[EdgeR]
Cufflinks
limma



  
  
 Visaualization of DEGs as Heatmaps 

Volcano Plots can be used to visalize significant chnges in gene expression




