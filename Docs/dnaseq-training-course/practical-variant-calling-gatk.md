##Genome Analysis Toolkit (GATK)

GATK is a collection of command-line tools for anlayzing high-throughput sequeccing data with a primary focus on variant discovery.The GATK pipeline workflow was applied following [best practices](https://software.broadinstitute.org/gatk/best-practices).
Variant Discovery

###HaplotypeCaller and Mutect2

HaplotypeCaller is designed to call germline variants, while Mutect2 is designed to call somatic variants.


Germline variants are variants against the reference.

Somatic variants contrast between two samples against the reference.
```
gatk HaplotypeCaller \
    -R reference.fasta \
    -I sample1.bam \
    -O variants.g.vcf \
    -ERC GVCF
    
```


GVCF stands for Genomic VCF. A GVCF is a kind of VCF, so the basic format specification is the same as for a regular VCF (see the spec documentation here), but a Genomic VCF contains extra information.

##Step 2: Prepare Analysis Ready Reads

####1. Convert sam to bam
####2. Sort BAM
####3. Mark Duplicate Reads

Check statistics using samtools flagstat

###For Germline variant calling

```
#########################################################################
for sample in $DIR/NIST7035_TAAGGCGA_Merged_001_trimmed_aln_PE.sam

do

echo $sample  
   filename=$(echo ${sample} | sed 's/.sam//')  
   echo $filename  
   
#1. Convert file from SAM to BAM format  
   
samtools view -b $sample > ${filename}.bam  
   
#2. Sort BAM file  
   
samtools sort ${filename}.bam -o ${filename}.sorted.bam   
   
#Picard tool to correct the header
#java -Xms128m -Xmx1024m -jar /cm/shared/apps/picard/2.6.0/picard.jar AddOrReplaceReadGroups I=${filename}.sorted.bam O=${filename}.picard.sorted.bam SORT_ORDER=coordinate RGID=foo RGLB=bar RGPL=Illumina RGSM=${filename} RGPU=JinliLuo
   
#3.Picard tool to mark duplicates
java -Xms128m -Xmx1024m -jar /cm/shared/apps/picard/2.6.0/picard.jar MarkDuplicates \
      I=${filename}.picard.sorted.bam \
      O=${filename}.markdup.picard.sorted.bam \
      M=${filename}.marked_dup_metrics.txt   

#4.index the bam file  
  
samtools index ${filename}.markdup.picard.sorted.bam  

#5.Check statistics using samtools flagstat

samtools flagstat ${filename}.markdup.picard.sorted.bam  > ${filename}.markdup.picard.sorted.bam.stat

   
#6. Single-sample GVCF calling with allele-specific annotations

gatk --java-options "-Xmx4g" HaplotypeCaller -R /data/crukjl19/reference/human_g1k_v37.fasta -I  ${filename}.markdup.picard.sorted.bam  -ERC GVCF -O ${filename}.markdup.picard.sorted.gatk.raw.g.vcf  -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation


done

################################################################################################
```
###Using samtools flagstat command to check the statistics of sorted BAM file
```
samtools flagstat NIST7035_TAAGGCGA_Merged_001_trimmed_aln_PE.sorted.bam
```

```
39800852 + 0 in total (QC-passed reads + QC-failed reads)
39789704 + 0 primary
0 + 0 secondary
11148 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
39774241 + 0 mapped (99.93% : N/A)
39763093 + 0 primary mapped (99.93% : N/A)
39789704 + 0 paired in sequencing
19894852 + 0 read1
19894852 + 0 read2
39527266 + 0 properly paired (99.34% : N/A)
39750284 + 0 with itself and mate mapped
12809 + 0 singletons (0.03% : N/A)
28872 + 0 with mate mapped to a different chr
13550 + 0 with mate mapped to a different chr (mapQ>=5)
```
###Mark duplicate reads

The aim of this step is to locate and tag duplicate reads in the BAM file.  
Duplicate reads are defined as originating from a a single fragment of DNA.  It can arise during smaple preparattion.

```
39800852 + 0 in total (QC-passed reads + QC-failed reads)
39789704 + 0 primary
0 + 0 secondary
11148 + 0 supplementary
2545449 + 0 duplicates
2545449 + 0 primary duplicates
39774241 + 0 mapped (99.93% : N/A)
39763093 + 0 primary mapped (99.93% : N/A)
39789704 + 0 paired in sequencing
19894852 + 0 read1
19894852 + 0 read2
39527266 + 0 properly paired (99.34% : N/A)
39750284 + 0 with itself and mate mapped
12809 + 0 singletons (0.03% : N/A)
28872 + 0 with mate mapped to a different chr
13550 + 0 with mate mapped to a different chr (mapQ>=5)
```

There are 2545449 are duplicate reads


3. Base Quality Score Recalibration (BQSR)

The last step of pre-processing mapped reads is the base quality score recalibration (BQSR) stage. The GATK tools detects systematic errors made by the sequencing machine while estimating the accuracy of each base. The systematic errors can have various sources ranging from technical machine errors to the variability in the sequencing chemical reactions. The two step BQSR process applies machine learning to model the possible errors and adjust the base quality scores accordingly. More details here. #

Analysis
The sequencing data were aligned by bwa mem6 against b37 human decoy reference genome. The alignments were sorted and PCR duplicates were marked by Picard (http://picard.sourceforge.net). For AJ trios, a joint variant calling was performed by GATK7 HaplotypeCaller on all three samples. For the Chinese son, both single sample variant calling (a VCF file) and the first step in cohort analysis (a gVCF file) were performed by GATK HaplotypeCaller. All variants in VCF files were quality filtered by standard GATK SNP variant quality score recalibration and indel hard filtration according to GATK Best Practices recommendations8,9.





Somatic short variant discovery (SNVs +Indels)

Identify somatic short variants (SNVs and Indels) in one or more tumor samples from a single individual, with or without a matched normal sample

###Expeacted Input

Bam files for each input tumor and normal samples

###Call candidate variants


> - [dna-alignment-snpcalling-bwa-gatk.md](dna-alignment-snpcalling-bwa-gatk.md) will be generated in the folder after the run is completed

