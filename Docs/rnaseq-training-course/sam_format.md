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
SAM files are TSV (tab-separated-values) files and begin with an optional header. The header consists of multiple lines, starting with an '@' character, each line is a record. 

Each record starts with its identifier and is followed by tab-separated tags. 

Each tag in the header consists of a two-character identifier, followed by ':', followed by the value.

|Record Name	|Tag 1 |Description | Tag2 | Description|
| ----------- |---------------|----------- |---------------|--------|
|@HD	|VN|	SAM version|SO	|sort order|
|@SQ	|SN|	the reference sequence names| LN| sequence lengths|

There also are other header record types.

The optional header section is followed by the alignment records. The alignment records are again tab-separated. There are 11 mandatory columns.


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


BAM files store their header as plain-text SAM headers. However, they additionally store the name and length information about the reference sequences. This information is mandatory since in BAM, the alignment records only contain the numeric ids of the reference sequences. Thus, the name is stored outside the record in the header.

> :memo: **title**

    
    > -  list under lists
    > -  under lists

        


> :bulb: if your editor gets confused by
not having and enclosing * then
just add it to end of abbr def.

---

> :warning: Don't indent these, doesn't seem to work

> :memo: **Memo Admonition**
        
use blockquotes
with emoji indicators for
admonition memos, callout etc..

---

> :boom:
> 
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
