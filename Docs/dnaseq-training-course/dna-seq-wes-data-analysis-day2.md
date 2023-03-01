#DNA-seq (WES) Data Analysis


### [Day2](dna-seq-wes-data-analysis-day2.md)

* Introduction of Genetic Variants
* Identify simple variants using GATK HaplotypeCaller 
* Visualise simple variant data (VCF files) 
* Perform basic variant filtering



###Introduction of Genetic Variation

There are approximately 3 billion base pairs in the human genome. Humans share 99.9% of DNA with other humans.  It is around 0.1% which is responsible for all the differenes that make each one of us individual.

Genetic variation is the difference in DNA sequences between individuals within a population. A gene variant is a permanent change in the DNA sequence that makes up a gene. 

> Single nucleotide polymorphisms (SNP)

The most common form of genetic variants among individuals are the smallest, known as single nucleotide polymorphisms (SNPs), describing a change in a single nucleotide anywhere in the genome. A nucleotide change is considered an SNP if the changes at a particular position are seen in more than 1% of the population.

As SNPs occur on average every 300 nucleotides, there will be appromimately 10 million SNPs in a person's genome.

The majority of SNP are harmless.  Some SNPs are the direct cause of a particular condition.

> Classification of genetic variants based on cell type or inheritance:  

![Based on Cell Type:](/assets/germline_somatic_mutation.jpg)

Germline Variants: are Inherited (or hereditary) variants are passed from parent to child and are present throughout a person’s life in virtually every cell in the body. These variants are also called germline variants because they are present in the parent’s egg or sperm cells, which are also called germ cells. 

Somatic Variants: are Non-inherited variants occur at some time during a person’s life and are present only in certain cells, not in every cell in the body. Because non-inherited variants typically occur in somatic cells (cells other than sperm and egg cells), they are often referred to as somatic variants. These variants cannot be passed to the next generation. Non-inherited variants can be caused by environmental factors such as ultraviolet radiation from the sun or can occur if an error is made as DNA copies itself during cell division.  Some somatic variants will lead to the development of diseases such as cancer. 

Both germline and somatic variants have been implicated in many cancers. 


>Classification by Effect on Sequence and Structure 

Ensembl Sequence Variants (Small-scale mutations)

|Type	|Description	|Example (Reference / Alternative)|
| ----------- | --------------|-------------------------------------------------|
|SNP	|Single Nucleotide Polymorphism	|Ref: ...TTG<span style="color:blue">**A**</span>CGTA...	
|||Alt: ...TTG<span style="color:red">**G**</span>CGTA...|
|Insertion	|Insertion of one or several nucleotides	|Ref: ...TTG<span style="color:blue">**AC**</span>GTA...	|
|||Alt: ...TTGA<span style="color:red">**TG**</span>CGTA...|
|Deletion	|Deletion of one or several nucleotides	|Ref: ...TTG<span style="color:blue">**A**</span>CGTA...|
|||Alt: ...TTGGTA...|
|Indel	|An insertion and a deletion, affecting 2 or more nucleotides	|Ref: ...TTG<span style="color:blue">**A**</span>CGTA...	|
|||Alt: ...TTG<span style="color:red">**GCT**</span>CGTA...|
|Substitution	|A sequence alteration where the length of the change in the variant is the same as that of the reference.	|Ref: ...TTG<span style="color:blue">**AC**</span>GTA...	|
|||Alt: ...TTG<span style="color:red">**TA**</span>GTA...|

Structural Variants (Large-scale mutations)

|Type	|Description	|Example (Reference / Alternative)|
| ----------- | --------------|-------------------------------------------------|
|CNV	|Copy Number Variation: increases or decreases the copy number of a region of DNA|Reference:"Gain" of one copy:![gene-duplication](/assets/gene-duplication.jpg)|
|||"Loss" of one copy:![gene-deletion](/assets/deletion_of_chromosome_section.jpg)|
|Inversion	|A continuous nucleotide sequence is inverted in the same position	|Reference:![gene-inversion](/assets/inversion_mutation.jpg)
Alternative:|
|Translocation	|A region of nucleotide sequence that has translocated to a new position|Reference:![gene-translocation](/assets/translocation_mutation2.jpg)
Alternative:|



Copy number variants can result in too much or too little of a protein being produced. Tumour cells often have structural variants that result in copy number variants. Individuals with other disorders such as epilepsy can also have copy number variants of genes associated with the disorder.





#### Mutations and recombination are major sources of variation.

Mutations are the original source of genetic variation. A mutation is a permanent alteration to a DNA sequence. De novo (new) mutations occur when there is an error during DNA replication that is not corrected by DNA repair enzymes. It is only once the error is copied by DNA replication, and fixed in the DNA that it is considered to be a mutation ![mutation](/assets/mutation2.jpg). Mutations may be beneficial to the organism; deleterious (harmful) to the organism; or neutral (have no effect on the fitness of the organism). 

Recombination is another major source of genetic variation Each of us has a mixture of genetic material from our parents. The mixing of this genetic material occurs during recombination when homologous DNA strands align and cross over. Recombination effectively ‘shuffles’ maternal and paternal DNA, creating new combinations of variants in the daughter germ-cells (Figure).
![Homologous Recombination ](/assets/homologous_recombination.jpg)

####WGS mainly focus on two types of variants

Single nucleotide polymorphisms (SNPs) and Short insertions and eletions (indels).

A SNP involves variation in just one nucleotide in a DNA sequence.  SNPs account for most of our genetic differences.Different populations may have different allele frequencies for polymorphic genes. Other types of variations such as deletions and insertions of nucleotides in DNA sequences account for a much smaller proportion of our overall gentic variation.

Tools available: GATK, Strelka, SppedSeq, Samtools, Varscan2, DRAGEN and DeepVariant.


####Variant identification and analysis

> Variant Calling

Carry out whole genome or whole exome sequencing to create FASTQ files.
Align the sequences to a reference genome, creating BAM or CRAM files.
Identify where the aligned reads differ from the reference genome and write to a VCF file.

> Somatic versus germline variant calling

In germline variant calling, the reference genome is the standard for the species of interest. This allows us to identify genotypes. As most genomes are diploid, we expect to see that at any given locus, either all reads have the same base, indicating homozygosity, or approximately half of all reads have one base and half have another, indicating heterozygosity. An exception to this would be the sex chromosomes in male mammals.

In somatic variant calling, the reference is a related tissue from the same individual. Here, we expect to see mosaicism between cells.

> Understanding VCF format

VCF is the standard file format for storing variation data. It is used by large scale variant mapping projects such as The International Genome Sample Resource(IGSR) https://www.internationalgenome.org/). It is also the standard output of variant calling software such as GATK and the standard input for variant analysis tools such as the Ensembl Variant Effect Predictor (VEP) (http://www.ensembl.org/info/docs/tools/vep/index.html) or for variation archives like European Variation archive (EVA) (https://www.ebi.ac.uk/eva/).

VCF is a preferred format because it is unambiguous, scalable and flexible, allowing extra information to be added to the info field. Many millions of variants can be stored in a single VCF file. 

> Variant Identifiers

Variants may have identifiers from multiple databases. You will see these different types of identifiers used throughout the literature and in other databases. Different types of identifiers are used for short variants and structural variants. Some common databases and examples of the identifiers they use are shown in the tables below.



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



### 
[Practical 3: Variant Calling](practical-variant-calling-gatk.md)


### 

