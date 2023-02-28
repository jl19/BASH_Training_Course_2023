##Practical:

###Step 1 - Copy and check the read data files

Example data used here is from Illumina genomic data of Pseudomonas aeruginosa (616 MB) (https://digitalinsights.qiagen.com/downloads/example-data/). The data is available from the
Short Read Archive under the Study Accession SRP010152.

The paired end data files can be found on SPECTRE at the following location:

/data/bioinf/Teaching/2019_NGS_Course/Data_QC/paeruginosa-reads/

Paired end 1 data file	= SRR396636.sra_1.fastq
Paired end 2 data file	= SRR396636.sra_2.fastq

![Screenshot of NoMachine login screen](/assets/NoMachine_Scratchfolder.png)


That paired end data is contained in two separate files, often known as Left and Right reads.
The first task is to copy the fastq data files to your scratch directory so you can use them.  Use the ‘cd’ (change directory) command to move to your scratch directory.  ‘cd $SCRATCHDIR’ takes you directly to your scratch directory, instead of having to type the full path which is in this format: /scratch/a/ab123.  Another advantage of $SCRATCHDIR is that you can just type ‘cd $SC’ and hit the tab key which will auto-complete the full command

Use the command below, followed by pressing the tab key to auto-complete, then hit return to move to your scratch directory

```
$  cd $SC 

```

The command ‘pwd’ which stands for ‘print working directory’ will tell you the path of your current directory. 
```
$  pwd
```
The output from the ‘pwd’ command will look something like this, (except your own 
username will come after the /scratch/ part)
/scratch/bbash/jl19



In order to keep things tidy in your scratch directory it is good practise to create separate directories for the different exercises.  We will now create a ‘Data_QC’ directory in your scratch directory.  The ‘mkdir’ command, short for ‘make directory’ will create a directory with the name you specify.

```
$  mkdir Data_QC
```

Then ‘cd’ to the Data_QC directory you just created.

```
$  cd Data_QC
```

We will now use the ‘cp’ command to copy the fastq data files from their current location (/data/bioinf/Teaching/2019_NGS_Course/Data_QC/) to your location (the Data_QC directory you just created).  Once you have run the ‘cp’ command use the ‘ls’ (list) command to list the contents of your directory and check that the files have been copied to your directory.  


```
$  cp /data/bioinf/Teaching/2017_NGS_Course/Data_QC/paeruginosa-reads/*  .
  
$  ls
```
NOTE: The ‘*’ at the end of the command tells ‘cp’ to copy everything in the Data_QC directory.  The ‘.’ tells the ‘cp’ command to copy the fastq files to the current location i.e. your Data_QC directory.

You will notice that the ‘cp’ command copied 2 fastq files and another file called ‘TruSeq3-PE-2.fa’ this is a file containing adapter sequences that we will use later to remove adapter contamination from the read data.

Use the ‘more’ command to view one of your fastq files
```
$  more SRR396637.sra_1.fastq
```
To exit the ‘more’ command and return to the command prompt press ‘q’


Try the ‘head’ command which displays just the first 10 lines of a file
```
$  head SRR396637.sra_1.fastq

```

Checked the fastq files are copied correctly. 

###Step 2: Assess the quality of the data using FastQC

FastQC is written by Simon Andrews of Babraham Bioinformatics, is a very popular quality control tool used to provide an overview of basic quality control metrics for raw next generation sequencing data.


The FastQC software needs to be made accessible to you on SPECTRE by using the ‘module load’ command

```
$  module load fastqc/0.11.5
```
Then launch FastQC as follows:
```
$  fastqc &
```
![Screenshot of fastqc](/assets/NoMachine_fastqc.png)

Run FastQc on command line, [SRR396636.sra_1_fastqc.html](/assets/SRR396636.sra_1_fastqc.html) will be generated in the folder after the run is completed

```
[jl19@spectre10 paeruginosa-reads]$ fastqc SRR396636.sra_1.fastq 
Picked up JAVA_TOOL_OPTIONS: -XX:MaxHeapSize=2048m
Started analysis of SRR396636.sra_1.fastq
Approx 5% complete for SRR396636.sra_1.fastq
Approx 10% complete for SRR396636.sra_1.fastq
Approx 15% complete for SRR396636.sra_1.fastq
Approx 20% complete for SRR396636.sra_1.fastq
Approx 25% complete for SRR396636.sra_1.fastq
Approx 30% complete for SRR396636.sra_1.fastq
Approx 35% complete for SRR396636.sra_1.fastq
Approx 40% complete for SRR396636.sra_1.fastq
Approx 45% complete for SRR396636.sra_1.fastq
Approx 50% complete for SRR396636.sra_1.fastq
Approx 55% complete for SRR396636.sra_1.fastq
Approx 60% complete for SRR396636.sra_1.fastq
Approx 65% complete for SRR396636.sra_1.fastq
Approx 70% complete for SRR396636.sra_1.fastq
Approx 75% complete for SRR396636.sra_1.fastq
Approx 80% complete for SRR396636.sra_1.fastq
Approx 85% complete for SRR396636.sra_1.fastq
Approx 90% complete for SRR396636.sra_1.fastq
Approx 95% complete for SRR396636.sra_1.fastq
Analysis complete for SRR396636.sra_1.fastq
[jl19@spectre10 paeruginosa-reads]$ 

```

Interpreation of the QC Report

FastQC has a really well documented [manual page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) with detailed explanations about every plot in the report

Trim the reads
```
[jl19@spectre10 paeruginosa-reads]$ module load skewer/0.2.2

```

```
[jl19@spectre10 paeruginosa-reads]$ skewer -m pe -q 25 -Q 20 -n -z SRR396636.sra_1.fastq SRR396636.sra_2.fastq
.--. .-.
: .--': :.-.
`. `. : `'.' .--. .-..-..-. .--. .--.
_`, :: . `.' '_.': `; `; :' '_.': ..'
`.__.':_;:_;`.__.'`.__.__.'`.__.':_;
skewer v0.2.2 [April 4, 2016]
Parameters used:
-- 3' end adapter sequence (-x):	AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
-- paired 3' end adapter sequence (-y):	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
-- maximum error ratio allowed (-r):	0.100
-- maximum indel error ratio allowed (-d):	0.030
-- mean quality threshold (-Q):		20
-- end quality threshold (-q):		25
-- minimum read length allowed after trimming (-l):	18
-- file format (-f):		Sanger/Illumina 1.8+ FASTQ (auto detected)
Wed Jul 27 11:29:25 2022 >> started
|=================================================>| (100.00%)
Wed Jul 27 11:30:23 2022 >> done (57.594s)
1909263 read pairs processed; of these:
    125 ( 0.01%) degenerative read pairs filtered out
  86590 ( 4.54%) read pairs filtered out by quality control
   8684 ( 0.45%) short read pairs filtered out after trimming by size control
   3823 ( 0.20%) empty read pairs filtered out after trimming by size control
1810041 (94.80%) read pairs available; of these:
1643835 (90.82%) trimmed read pairs available after processing
 166206 ( 9.18%) untrimmed read pairs available after processing
log has been saved to "SRR396636.sra_1-trimmed.log".
[jl19@spectre10 paeruginosa-reads]$ 

```

Run FastQc on command line after trimming, [SRR396636.sra_1-trimmed-pair1_fastqc.html](/assets/SRR396636.sra_1-trimmed-pair1_fastqc.html) will be generated in the folder after the run is completed

Sequence Length Distribution

In many cases this will produce a simple graph showing a peak only at one size, but for variable length FastQ files this will show the relative amounts of each different size of sequence fragment.

Warning
This module will raise a warning if all sequences are not the same length.

Per Tile Sequence Quality


Reasons for seeing warnings or errors on this plot could be transient problems such as bubbles going through the flowcell, or they could be more permanent problems such as smudges on the flowcell or debris inside the flowcell lane.

Warning
This module will issue a warning if any tile shows a mean Phred score more than 2 less than the mean for that base across all tiles.

Here is a introduction on how to interpret FastQC report from Babraham Bioiformatics
[Using FastQC to check the quality of high throughput sequence](https://www.youtube.com/watch?v=bz93ReOv87Y)
