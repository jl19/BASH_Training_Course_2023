DNA variant calling

Genetic Variation

Genetic variation is the difference in DNA sequences between individuals within a population. Variation occurs in germ cells i.e. sperm and egg, and also in somatic (all other) cells. Only variation that arises in germ cells can be inherited from one individual to another and so affect population dynamics, and ultimately evolution. Mutations and recombination are major sources of variation.

Mutations are the original source of genetic variation. A mutation is a permanent alteration to a DNA sequence. De novo (new) mutations occur when there is an error during DNA replication that is not corrected by DNA repair enzymes. It is only once the error is copied by DNA replication, and fixed in the DNA that it is considered to be a mutation (Figure 1). Mutations may be beneficial to the organism; deleterious (harmful) to the organism; or neutral (have no effect on the fitness of the organism). 

Recombination is another major source of genetic variation Each of us has a mixture of genetic material from our parents. The mixing of this genetic material occurs during recombination when homologous DNA strands align and cross over. Recombination effectively ‘shuffles’ maternal and paternal DNA, creating new combinations of variants in the daughter germ-cells (Figure).
![Homologous Recombination ](/assets/homologous_recombination.jpg)

WGS mainly focus on two types of variants

Single nucleotide polymorphisms (SNPs) and Short insertions and eletions (indels).

Tools available: GATK, Strelka, SppedSeq, Samtools, Varscan2, DRAGEN and DeepVariant.


Variant identification and analysis

Variant Calling

Carry out whole genome or whole exome sequencing to create FASTQ files.
Align the sequences to a reference genome, creating BAM or CRAM files.
Identify where the aligned reads differ from the reference genome and write to a VCF file.

Somatic versus germline variant calling
In germline variant calling, the reference genome is the standard for the species of interest. This allows us to identify genotypes. As most genomes are diploid, we expect to see that at any given locus, either all reads have the same base, indicating homozygosity, or approximately half of all reads have one base and half have another, indicating heterozygosity. An exception to this would be the sex chromosomes in male mammals.

In somatic variant calling, the reference is a related tissue from the same individual. Here, we expect to see mosaicism between cells.

Understanding VCF format
VCF is the standard file format for storing variation data. It is used by large scale variant mapping projects such as The International Genome Sample Resource(IGSR) https://www.internationalgenome.org/). It is also the standard output of variant calling software such as GATK and the standard input for variant analysis tools such as the Ensembl Variant Effect Predictor (VEP) (http://www.ensembl.org/info/docs/tools/vep/index.html) or for variation archives like European Variation archive (EVA) (https://www.ebi.ac.uk/eva/).

VCF is a preferred format because it is unambiguous, scalable and flexible, allowing extra information to be added to the info field. Many millions of variants can be stored in a single VCF file. 

Variant Identifiers

Variants may have identifiers from multiple databases. You will see these different types of identifiers used throughout the literature and in other databases. Different types of identifiers are used for short variants and structural variants. Some common databases and examples of the identifiers they use are shown in the tables below.

| File Name   |File Extension | File Type  | Description                                      |                       
| ----------- | --------------|------------|--------------------------------------------------|
| `Fasta`     | .fasta, .fa   |sequences   |txt file for nucleotie or peptie sequencese       |
| `Fastq`     | .fastq, .fq   |read data   |txt file storing both sequence and its quality scores |
| `SAM`       | .sam          |short read alignment|sequence alignment file                   |
| `BAM`       |.bam           | binary SAM |
| `VCF`       | .vcf          | Variant information, Variant Call Format| txt file for variants calls incluing aNP, indels, CNV
|`GFF/GTF`    | .gff, .gtf    | annotation data | General Transfer/Feature Format annotation file



Table 1 Types of variant identifiers

|Identifier type	|Example|	Description|
| ----------- | --------------|-------------------------------------------------|
 |ssID	|ss335	|Submitted SNP ID assigned by dbSNP or EVA.|
 |rsID	|rs334	|Reference SNP ID assigned by dbSNP or EVA. ssIDs of the same variant type that colocalise are combined to give an rsID for that locus.|
| HGVS*	|ENST00000366667.4:c.803T>C	|Expresses the location of the variant in terms of a transcript or protein.|
| COSMIC ID	|COSM1290	|ID assigned by COSMIC for somatic variants.|
| HGMD	|CD830010	|ID assigned by HGMD to variants known to be associated with human inherited diseases.|
| ClinVar	|RCV000016573	|ID assigned to dbSNP or dbVar/DGVa annotated variants, linking them to human health.|
| UniProt	|VAR_010085	|ID assigned by UniProt for reviewed human.|
|DGVa variant call	|essv8691751	|Submitted structural variant ID assigned by DGVa. Variants are shared with dbVar.|
|dbVar variant call	|nssv1602417	|Submitted structural variant ID assigned by dbVar. Variants are shared with DGVa.|
|DGVa variant region	|esv3364878	|Variant region variant ID assigned by DGVa. Overlapping submitted variants (essv and nssv) are combined into a single variant region. The boundaries of a variant region may not match those of the submitted variants, which can vary.|
|dbVar variant region	|nsv916030	|Variant region variant ID assigned by dbVar. Overlapping submitted variants (essv and nssv) are combined into a single variant region. The boundaries of a variant region may not match those of the submitted variants, which can vary.|




VCF files are tab delimited text files. Here is an example of a variant in VCF (Figure 12) as viewed in a spreadsheet:
