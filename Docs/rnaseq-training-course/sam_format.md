* SAM / BAM File Structure
We will give an quick overview of the SAM and BAM formats here. Note that this overview serves more to remind you what the formats are about and are not meant to teach how to use the SAM and BAM format.

The following shows an example of a SAM file.
```
@HD VN:1.3  SO:coordinate
@SQ SN:ref  LN:45
@SQ SN:ref2 LN:40
r001    163 ref 7   30  8M4I4M1D3M  =   37  39  TTAGATAAAGAGGATACTG *   XX:B:S,12561,2,20,112
r002    0   ref 9   30  1S2I6M1P1I1P1I4M2I  *   0   0   AAAAGATAAGGGATAAA   *
r003    0   ref 9   30  5H6M    *   0   0   AGCTAA  *
r004    0   ref 16  30  6M14N1I5M   *   0   0   ATAGCTCTCAGC    *
r003    16  ref 29  30  6H5M    *   0   0   TAGGC   *
r001    83  ref 37  30  9M  =   7   -39 CAGCGCCAT   *
  
```
  SAM files are TSV (tab-separated-values) files and begin with an optional header. The header consists of multiple lines, starting with an '@' character, each line is a record. Each record starts with its identifier and is followed by tab-separated tags. Each tag in the header consists of a two-character identifier, followed by ':', followed by the value.

If present, the @HD record must be the first record and specifies the SAM version (tag VN) used in this file and the sort order (SO). The optional @SQ header records give the reference sequence names (tag SN) and lengths (tag LN). There also are other header record types.

The optional header section is followed by the alignment records. The alignment records are again tab-separated. There are 11 mandatory columns.

|Phred Quality Score          |Probability of Incorrect Base Call |   Base Call Accuracy|
| ----------- |---------------|--------------------------------------------------|
|10                            |  1 in 10                       |      90%|
|20                            |  1 in 100                       |     99%|
|30                            |  1 in 1,000                     |     99.9%|
|40                            |  1 in 10,000                     |    99.99%|
|50                            |  1 in 100,000                    |    99.999%|




|Col|	Field	Type|	N/A Value	|Description|
| ----------- |---------------|------------------------|--------------------------|
|1	|QNAME|	string|	mandatory	The query/read name.|
|2	|FLAG	|int|	mandatory	The record’s flag.|
|3	|RNAME|	string|	*	|The reference name.|
|4	|POS	|32-bit int|	0	|1-based position on the reference.|
|5	|MAPQ	|8-bit int|	255	|The mapping quality.|
|6	|CIGAR	|string	|*	|The CIGAR string of the alignment.|
|7|	RNEXT|	string|	*	|The reference of the next mate/segment.|
 |8	|PNEXT|	string|	0	|The position of the next mate/seqgment.|
|9|	TLEN|	string|	0	|The observed length of the template.|
|10	|SEQ|	string|	*|	The query/read sequence.|
|11|	QUAL|	string|	*	|The ASCII PHRED-encoded base qualities.|



Notes:
  
*  The SAM standard talks about “queries”. In the context of read mapping, where the format originates, queries are reads.
* The SAM standard talks about “templates” and “segments”. In the case of paired-end and mate-pair mapping the template consists of two segments, each is one read. The template length is the insert size.
* Paired-end reads are stored as two alignments records with the same QNAME. The first and second mate are discriminated by the FLAG values.
* When the FLAG indicates that SEQ is reverse-complemented, then QUAL is reversed.
* Positions in the SAM file are 1-based. When read into a BamAlignmentRecord (see below), the positions become 0-based.
* The qualities must be stored as ASCII PRED-encoded qualities.
* The query and reference names must not contain whitespace. It is common to trim query and reference ids at the first space.
* There are many ambiguities, recommendations, and some special cases in the formats that we do not describe here. We recommend that you follow this tutorial, start working with the SAM and BAM formats and later read the SAM specification “on demand” when you need it.

The 11 mandatory columns are followed by an arbitrary number of optional tags. Tags have a two-character identifier followed by ":${TYPE}:", followed by the tag’s value.

BAM files store their header as plain-text SAM headers. However, they additionally store the name and length information about the reference sequences. This information is mandatory since in BAM, the alignment records only contain the numeric ids of the reference sequences. Thus, the name is stored outside the record in the header.

A First Working Example
The following example shows an example of a program that reads the file with the path example.sam and prints its contents back to the user on stdout. If you want to try out this program then create a file with the sample SAM content from above and adjust the path "example.sam" in the program below to the path to your SAM file (e.g. "path/to/my_example.sam").


> :memo: **title**
    >
    > - list under lists
    > - under lists

> :memo: **title**
    >
    > - list under lists
    > - under lists


> :bulb: if your editor gets confused by
not having and enclosing * then
just add it to end of abbr def.

---

>:warning: Don't indent these, doesn't seem to work

> :memo: **Memo Admonition**
use blockquotes
with emoji indicators for
admonition memos, callout etc..

---

> :boom:
Title title like above is optional

---

> :bulb: See [the section about blocks](blocks.md#cheatsheet)
for the list of emojis that can be used.
Memo Admonition

use blockquotes with emoji indicators for admonition memos, callout etc..

Title title like above is optional

See the section about blocks for the list of emojis that can be used.






  # tutorial_basic_sam_bam_io_example1
  @HD     VN:1.3  SO:coordinate
@SQ     SN:ref  LN:45
@SQ     SN:ref2 LN:40
r001    163     ref     7       30      8M4I4M1D3M      =       37      39      TTAGATAAAGAGGATACTG     *       XX:B:S,12561,2,20,112
r002    0       ref     9       30      1S2I6M1P1I1P1I4M2I      *       0       0       AAAAGATAAGGGATAAA       *
  r003    0       ref     9       30      5H6M    *       0       0       AGCTAA  *
  r004    0       ref     16      30      6M14N1I5M       *       0       0       ATAGCTCTCAGC    *
  r003    16      ref     29      30      6H5M    *       0       0       TAGGC   *
  r001    83      ref     37      30      9M      =       7       -39     CAGCGCCAT       *
  83      *       *       *       *       *       0       *       *       *
  83      *       *       *       *       *       0       *       *       *
  ...