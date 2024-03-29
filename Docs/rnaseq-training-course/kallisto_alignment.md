# Pseudoalignment Using kallisto

kallisto is a program for quantifying abundances of transcripts from bulk and single-cell RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads.
It is based on the novel idea of pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment.


kallisto is fast and also can be used as quantification tool.

## Download the cDNA reference and build a index file.

### Download the cDNA reference
```
#Go to the directory

cd /scratch/bbash/jl19/Data_QC/mm10/

#Command wget stands for web get, use it to download cDNA ref from a ftp server

wget ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz
```
### Load the kallisto module
```
module load kallisto/0.46.1
```

### Kallisto Index

```
kallisto index -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx '/scratch/bbash/jl19/Data_QC/mm10/Mus_musculus.GRCm38.cdna.all.fa.gz' 

#If already in the directory of /scratch/bbash/jl19/Data_QC/mm10/, do not need the path of directory

[jl19@spectre13 mm10]$ kallisto index -i mm_10_genome.idx Mus_musculus.GRCm38.cdna.all.fa.gz 

[build] loading fasta file Mus_musculus.GRCm38.cdna.all.fa.gz
[build] k-mer length: 31
[build] warning: clipped off poly-A tail (longer than 10)
        from 650 target sequences
[build] warning: replaced 3 non-ACGUT characters in the input sequence
        with pseudorandom nucleotides
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done 
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 738750 contigs and contains 101004965 k-mers 
```
Time Required: 3min40sec

### Pseudoaligment using kallisto quantification
```
cd /scratch/bbash/jl19/Data_QC/
```
or cd .. if in already in the directory  of /scratch/bbash/jl19/Data_QC/mm10/
```
cd ..
```
```
mkdir kallisto_output 

cd kallisto_output
```
```
# Replace /bbash/jl19/ with your own directory name
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457551_24_hours-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457551_24_hours-trimmed.fastq.gz;
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457552_24_hours-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457552_24_hours-trimmed.fastq.gz;
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457561_24_hours-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457561_24_hours-trimmed.fastq.gz;
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457562_24_hours-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457562_24_hours-trimmed.fastq.gz;

kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457553_2_hours-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457553_2_hours-trimmed.fastq.gz;
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457554_2_hours-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457554_2_hours-trimmed.fastq.gz;
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457555_2_hours-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457555_2_hours-trimmed.fastq.gz;
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457556_2_hours-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457556_2_hours-trimmed.fastq.gz;

kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457557_control-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457557_control-trimmed.fastq.gz;
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457558_control-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457558_control-trimmed.fastq.gz;
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457559_control-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457559_control-trimmed.fastq.gz;
kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457560_control-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457560_control-trimmed.fastq.gz;
```

```
[jl19@spectre13 kallisto_output]$ kallisto quant -i /scratch/bbash/jl19/Data_QC/mm10/mm_10_genome.idx -o /scratch/bbash/jl19/Data_QC/kallisto_output/SRR7457551_24_hours-trimmed -b 100 --single -l 100 -s 20  /scratch/bbash/jl19/Data_QC/trimmed/SRR7457551_24_hours-trimmed.fastq.gz;

[quant] fragment length distribution is truncated gaussian with mean = 100, sd = 20
[index] k-mer length: 31
[index] number of targets: 119,414
[index] number of k-mers: 101,004,965
[index] number of equivalence classes: 436,756
[quant] running in single-end mode
[quant] will process file 1: /scratch/bbash/jl19/Data_QC_May2023/trimmed/SRR7457551_24_hours-trimmed.fastq.gz
[quant] finding pseudoalignments for the reads ... done
[quant] processed 29,108,478 reads, 25,714,780 reads pseudoaligned
[   em] quantifying the abundances ... done
[   em] the Expectation-Maximization algorithm ran for 1,171 rounds
[bstrp] running EM for the bootstrap: 100
```
Time Required: 40 mins
output files from kallisto:

> - abundance.h5: a HDF5 binary file containing run info, abundance esimates, bootstrap estimates, and transcript length information length. This file can be read in by sleuth
> - abundance.tsv: is a plaintext file of the abundance estimates. It does not contains bootstrap estimates. Please use the --plaintext mode to output plaintext abundance estimates. Alternatively, kallisto h5dump can be used to output an HDF5 file to plaintext. The first line contains a header for each column, including estimated counts, TPM, effective length.
```
target_id            length  eff_length est_counts  tpm
ENSMUST00000178537.2    12      4.16231      0       0
ENSMUST00000178862.2    14      4.39563      0       0
ENSMUST00000196221.2    9       3.70903      0       0
ENSMUST00000179664.2    11      4.02702      0       0
ENSMUST00000177564.2    16      4.59059      0       0
ENSMUST00000179520.2    11      4.02702      0       0
ENSMUST00000179883.2    16      4.59059      0       0
ENSMUST00000040583.7    7763    7664         3542.85 29.3149
ENSMUST00000220369.2    4981    4882         38.7335 0.503131
ENSMUST00000218847.2    912     813          20.0819 1.56642
ENSMUST00000217966.2    780     681          2.60721 0.242785
ENSMUST00000219641.2    3657    3558         45.5654 0.812123
ENSMUST00000218186.2    2222    2123         22.737  0.679166
ENSMUST00000218254.2    4111    4012         50.9482 0.805305
ENSMUST00000219277.2    5052    4953         40.1902 0.514569
```


> - run_info.json: is a json file containing information about the run

All the output files can be copied from:
```
cp -r /data/bioinf/Teaching/2023_NGS_Course/Data_QC/RNA-Seq-GSE116583/kallisto_output/. ./
```
## [Back to Day 1 Practicals](rna-seq-wes-data-analysis-day1.md/#day-1-practicals)
