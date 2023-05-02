# RNA-Seq Data Pre-processing

>  Example data Used:
> Bioproject PRJNA450151
> SRA Study Number: SRP151689
> Title: Mus musculus Transcriptome or Gene expression
> Abstract:Alveolar macrophages were isolated in a murine model of lung transplant to study their role in ischemia reperfusion injury and primary graft dysfunction.

### Step 1 - Import/Copy and check the read data files

#### Importing raw data files

* Importing from [SRA Explorer](https://sra-explorer.info/)
* Using SRA study number (SRP151689) to search
* Select required datasets and add to collection
* Goes to saved datasets
* Datasets could be then downloaded using Raw FastQ Download URLs or Bash script for downloading FastQ files

*********
* Importing from [The European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/)
* Using SRA study number (SRP151689) search
* Select corresponding study
* Select required data/samples to download

>  [Offline HOW To download data can be found here](how-to-SRA-Explorer.md)

***********

#### Coping raw data files

The example raw data (single end data) files can also be copied from SPECTRE at the following location:

/data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/raw_data/

> - SRR7457557_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz
> - SRR7457558_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz
> - SRR7457559_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz
> - SRR7457560_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz
> - SRR7457553_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz
> - SRR7457554_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz
> - SRR7457555_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz
> - SRR7457556_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz
> - SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
> - SRR7457552_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
> - SRR7457561_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
> - SRR7457562_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz

![Screenshot of NoMachine login screen](https://jl19.github.io/BASH_Training_Course_2023/Docs/assets/NoMachine_Scratchfolder.png)

To copy the fastq data files to your scratch directory so you can use them.  Use the ‘cd’ (change directory) command to move to your scratch directory.  
```
cd $SCRATCHDIR
```
It takes you directly to your scratch directory, instead of having to type the full path which is in this format: /scratch/a/ab123.  
Type ‘cd $SC’ and press the tab key which will auto-complete the full command 'cd $SCRATCHDIR'
Use the command below, followed by pressing the tab key to auto-complete, then hit return to move to your scratch directory

```
cd $SC 
```
The command ‘pwd’ which stands for ‘print working directory’ will tell you the path of your current directory. 
```
pwd
```
The output from the ‘pwd’ command will look something like this, (except your own 
username will come after the /scratch/ part)
/scratch/bbash/jl19

In order to keep things tidy in your scratch directory it is good practise to create separate directories for the different exercises.  We will now create a ‘Data_QC’ directory in your scratch directory.  The ‘mkdir’ command, short for ‘make directory’ will create a directory with the name you specify.
```
mkdir Data_QC
```
Then ‘cd’ to the Data_QC directory you just created.

```
cd Data_QC
```
Use the ‘cp’ command to copy the fastq data files from their current location (/data/bioinf/Teaching/2022_NGS_Course/Data_QC/) to your location (the Data_QC directory you just created).   

Use the ‘ls’ (list) command to list the contents of your directory and check that the files have been copied to your directory.  

```
cp -r  /data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/raw_data/. ./
ls
```
Time Requied: 20s

> NOTE: The ‘.’ at the end of the command tells ‘cp’ to copy everything in the Data_QC directory.  The ‘.’ tells the ‘cp’ command to copy the fastq files to the current location i.e. your Data_QC directory.

The ‘cp’ command copied 12 fastq.gz files as follow:

```
[jl19@spectre12 Data_QC]$ls

SRR7457557_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz
SRR7457558_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz
SRR7457559_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz
SRR7457560_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz
SRR7457553_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz
SRR7457554_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz
SRR7457555_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz
SRR7457556_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz
SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
SRR7457552_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
SRR7457561_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
SRR7457562_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
```

Checked the fastq files are copied correctly. 

### Step 2: Assess the quality of the data using FastQC

FastQC is written by Simon Andrews of Babraham Bioinformatics, is a very popular quality control tool used to provide an overview of basic quality control metrics for raw next generation sequencing data.

To check software installed (available) on SPECTRE using the following command:
```
module ava
```
To load the specific software using the ‘module load’ command followed by the module name:

```
module load 'modulename'
```
The FastQC software can be accessed on SPECTRE by using the following command

```
module load fastqc/0.12.1
```
```
[jl19@spectre12 Data_QC]$ module load fastqc/0.12.1
Loading fastqc/0.12.1
Loading requirement: java/1.8
``` 

Run FastQc software:

```
fastqc &
```
![Screenshot of fastqc](https://jl19.github.io/BASH_Training_Course_2023/Docs/assets/fastqc.jpg)

>  [SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant_fastqc.html](https://jl19.github.io/BASH_Training_Course_2023/Docs/assets/SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant_fastqc.html) will be generated in the folder after the run is completed

Run FastQc on command line:
```  
[jl19@spectre12 Data_QC]$ fastqc SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz

Picked up JAVA_TOOL_OPTIONS: -XX:MaxHeapSize=2048m
[jl19@spectre13 Data_QC]$ fastqc SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
application/gzip
Started analysis of SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
Approx 5% complete for SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
Approx 10% complete for SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
Approx 15% complete for SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
Approx 20% complete for SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
Approx 25% complete for SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
....
Analysis complete for SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz
[jl19@spectre12 Data_QC]$ 


#output file: SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant_fastqc.zip
```
Time taken around 2.5 mins
> - FastQC has a really well documented [manual page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) with detailed explanations about every plot in the report

### Step 3. Trim the reads

```
#make diretory for trimmed files
mkdir trimmed

#load module
module load skewer/0.2.2

#list skewer USAGE
skewer --help
```
> - Output:
```
Skewer (A fast and accurate adapter trimmer for paired-end reads)
Version 0.2.2 (updated in April 4, 2016), Author: Hongshan Jiang

USAGE: skewer [options] <reads.fastq> [paired-reads.fastq]
    or skewer [options] - (for input from STDIN)

OPTIONS (ranges in brackets, defaults in parentheses):
 Adapter:
          -x <str> Adapter sequence/file (AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)
          -y <str> Adapter sequence/file for pair-end reads (AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA),
                   implied by -x if -x is the only one specified explicitly.
          -M, --matrix <str> File indicates valid adapter pairing (all-ones matrix).
          -j <str> Junction adapter sequence/file for Nextera Mate Pair reads (CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG)
          -m, --mode <str> trimming mode; 1) single-end -- head: 5' end; tail: 3' end; any: anywhere (tail)
                           2) paired-end -- pe: paired-end; mp: mate-pair; ap: amplicon (pe)
          -b, --barcode    Demultiplex reads according to adapters/primers (no)
 Tolerance:
          -r <num> Maximum allowed error rate (normalized #errors / length of aligned region) [0, 0.5], (0.1)
          -d <num> Maximum allowed indel error rate [0, r], (0.03)
                   reciprocal is used for -r, -e and -d when num > or = 2
          -k <int> Minimum overlap length for adapter detection [1, inf);
                   (max(1, int(4-10*r)) for single-end; (<junction length>/2) for mate-pair)
 Clipping:
          -c, --cut <int>,<int> Hard clip off the 5' leading bases as the barcodes in amplicon mode; (no)
          -e, --cut3            Hard clip off the 3' tailing bases if the read length is greater than
                                the maximum read length specified by -L; (no)
 Filtering:
          -q, --end-quality  <int> Trim 3' end until specified or higher quality reached; (0)
          -Q, --mean-quality <int> The lowest mean quality value allowed before trimming; (0)
          -l, --min <int> The minimum read length allowed after trimming; (18)
          -L, --max <int> The maximum read length allowed after trimming; (no limit)
          -n  Whether to filter out highly degenerative (many Ns) reads; (no)
          -u  Whether to filter out undetermined mate-pair reads; (no)
          -N, --fillNs Whether to replace trimmed bases with Ns (has no effect with 'b' or '-m mp'); (no)
 Input/Output:
          -f, --format <str>   Format of FASTQ quality value: sanger|solexa|auto; (auto)
          -o, --output <str>   Base name of output file; ('<reads>.trimmed')
          -z, --compress       Compress output in GZIP format (no)
          -1, --stdout         Redirect output to STDOUT, suppressing -b, -o, and -z options (no)
          --qiime              Prepare the "barcodes.fastq" and "mapping_file.txt" for processing with QIIME; (default: no)
          --quiet              No progress update (not quiet)
          -A, --masked-output  Write output file(s) for trimmed reads (trimmed bases converted to lower case) (no)
          -X, --excluded-output Write output file(s) for excluded reads (no)
 Miscellaneous:
          -i, --intelligent     For mate-pair mode, whether to redistribute reads based on junction information; (no)
          -t, --threads <int>   Number of concurrent threads [1, 32]; (1)

EXAMPLES:
          skewer -Q 9 -t 2 -x adapters.fa sample.fastq -o trimmed
          skewer -x AGATCGGAAGAGC -q 3 sample-pair1.fq.gz sample-pair2.fq.gz
          skewer -x TCGTATGCCGTCTTCTGCTTGT -l 16 -L 30 -d 0 srna.fastq
          skewer -m mp -i lmp-pair1.fastq lmp-pair2.fastq
          skewer -m ap --cut 0,6 --qiime -x forward-primers.fa -y reverse-primers.fa mix-pair1.fastq mix-pair2.fastq
```

> - Perform Trimming

```
skewer -m any -q 25 -Q 20 -n -z SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457551_24_hours
```

```
[jl19@spectre13 Data_QC_May2023]$ skewer -m any -q 25 -Q 20 -n -z SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457551_24_hours
.--. .-.
: .--': :.-.
`. `. : `'.' .--. .-..-..-. .--. .--.
_`, :: . `.' '_.': `; `; :' '_.': ..'
`.__.':_;:_;`.__.'`.__.__.'`.__.':_;
skewer v0.2.2 [April 4, 2016]
Parameters used:
-- 3' end adapter sequence (-x):	AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
-- maximum error ratio allowed (-r):	0.100
-- maximum indel error ratio allowed (-d):	0.030
-- mean quality threshold (-Q):		20
-- end quality threshold (-q):		25
-- minimum read length allowed after trimming (-l):	18
-- file format (-f):		Sanger/Illumina 1.8+ FASTQ (auto detected)
-- minimum overlap length for adapter detection (-k):	3
Fri Apr 28 11:47:57 2023 >> started
|=================================================>| (100.00%)
Fri Apr 28 11:58:57 2023 >> done (659.910s)
29763084 reads processed; of these:
    7663 ( 0.03%) degenerative reads filtered out
    4107 ( 0.01%) reads filtered out by quality control
       8 ( 0.00%) short reads filtered out after trimming by size control
  642828 ( 2.16%) empty reads filtered out after trimming by size control
29108478 (97.80%) reads available; of these:
 2790093 ( 9.59%) trimmed reads available after processing
26318385 (90.41%) untrimmed reads available after processing
log has been saved to "trimmed/SRR7457551_24_hours-trimmed.log".
[jl19@spectre13 Data_QC_May2023]$ 
```
Time Required:11 mins

Continue to trim other raw data files using bash script.

```
#trim_fastq_file.sh
#!/usr/bin/env bash

skewer -m any -q 25 -Q 20 -n -z  SRR7457557_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz -o trimmed/SRR7457557_control;
skewer -m any -q 25 -Q 20 -n -z  SRR7457558_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz -o trimmed/SRR7457558_control;
skewer -m any -q 25 -Q 20 -n -z  SRR7457559_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz -o trimmed/SRR7457559_control;
skewer -m any -q 25 -Q 20 -n -z  SRR7457560_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages.fastq.gz -o trimmed/SRR7457560_control;
skewer -m any -q 25 -Q 20 -n -z  SRR7457553_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457553_2_hours;
skewer -m any -q 25 -Q 20 -n -z  SRR7457554_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457554_2_hours;
skewer -m any -q 25 -Q 20 -n -z  SRR7457555_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457555_2_hours;
skewer -m any -q 25 -Q 20 -n -z  SRR7457556_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_2_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457556_2_hours;
skewer -m any -q 25 -Q 20 -n -z  SRR7457551_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457551_24_hours;
skewer -m any -q 25 -Q 20 -n -z  SRR7457552_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457552_24_hours;
skewer -m any -q 25 -Q 20 -n -z  SRR7457561_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457561_24_hours;
skewer -m any -q 25 -Q 20 -n -z  SRR7457562_RNA-Seq_of_mus_musculus_adult_male_alveolar_macrophages_24_hours_after_lung_transplant.fastq.gz -o trimmed/SRR7457562_24_hours;
```
All trimmed files can be copied from server.

cp /data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/raw_data/trimmed/. ./

> - Run FastQc on command line after trimming

> - [SRR7457551_24_hours.fastq-trimmed_fastqc.html](https://jl19.github.io/BASH_Training_Course_2023/Docs/assets/SRR7457551_24_hours-trimmed.fastqc.html) will be generated in the folder after the run is completed
-------------------

### Step 4. Understanding the FastQC Report
> - Sequence Length Distribution

In many cases this will produce a simple graph showing a peak only at one size, but for variable length FastQ files this will show the relative amounts of each different size of sequence fragment.

Warning: This module will raise a warning if all sequences are not the same length.

> - Per Tile Sequence Quality

Reasons for seeing warnings or errors on this plot could be transient problems such as bubbles going through the flowcell, or they could be more permanent problems such as smudges on the flowcell or debris inside the flow cell lane.

Warning: This module will issue a warning if any tile shows a mean Phred score more than 2 less than the mean for that base across all tiles.

> - Here is a introduction on how to interpret FastQC report from Babraham Bioiformatics
[Using FastQC to check the quality of high throughput sequence](https://www.youtube.com/watch?v=bz93ReOv87Y)

[Go to Day1 Practical 2: Align reads to a reference genome using HISAT2](align-reads-to-reference-genome.md)

## [Back to Day 1 Practicals](rna-seq-wes-data-analysis-day1.md#day-1-practicals)

