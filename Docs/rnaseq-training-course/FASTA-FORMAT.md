# FASTA Format

## FASTA format description (more information available from [NCBI(https://blast.ncbi.nlm.nih.gov/doc/blast-topics/)]

FASTA stands for "Fast All" or "FastA", which is a popular format for representing nucleotide or protein sequences and their associated metadata. The FASTA format was introduced by David J. Lipman and William R. Pearson in 1985 as a simple and efficient way to search sequence databases using sequence similarity.

FASTA format is a text-based format for representing either nucleotide sequences or protein sequences, in which base pairs or amino acids are represented using single-letter codes. A sequence in FASTA format begins with a single-line description, followed by lines of sequence data. The description line is distinguished from the sequence data by a greater-than (">") symbol in the first column. It is recommended that all lines of text be shorter than 80 characters in length.
An example sequence in FASTA format is:

FASTA format is characterized by its two-line structure. 

|Line|Description|
|-----|------|
|The first line |starts with a ">" symbol, followed by a unique identifier for the sequence, and then optional additional information about the sequence separated by space. |
|The second line |It contains the actual sequence of nucleotides or protein, with no spaces or line breaks.|

Here is an example of a FASTA-formatted sequence:
```
>gi|123456|ref|NM_001234.5| Homo sapiens mRNA for example gene, transcript variant 1, complete cds
ATGGCGACGCTGCGCGCGCTGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG

```
|Name|Description|
|---|----|
|gi| GenInfo Identifier|It is a unique numerical identifier assigned to a nucleotide or protein sequence record in the GenBank database. |

The GenBank database is a comprehensive collection of publicly available DNA and protein sequences maintained by the National Center for Biotechnology Information (NCBI).

When a protein sequence is obtained from the GenBank database, the sequence identifier in the FASTA format is usually prefixed with "gi|" followed by the numerical identifier of the record. For example, in the following FASTA format sequence identifier:

```
>gi|129361|sp|P01013|OVAX_CHICK GENE X PROTEIN (OVALBUMIN-RELATED)
```
"gi" indicates that the protein sequence is from the GenBank database, "129361" is the numerical identifier of the record, and the rest of the information provides details about the protein's origin, function, and other properties. However, note that the use of "gi" identifiers in FASTA format is less common now, as NCBI has phased out the use of "gi" identifiers in favor of accession numbers.

```
>sp|P69905|TBA1A_HUMAN Tubulin alpha-1A chain OS=Homo sapiens GN=TUBA1A PE=1 SV=2
MREIVHIQAGQCGNQIGAKFWEVISDEHGIDPSGNYVGDSDLQLERINVYYNEATGGKYVP
RAILVDLEPGTMDSVRSGPFGQIFRPDNFVFGQTGAGNNWAKGHYTEGAELVDSVLDVVR
KEAESCDCLQGFQLTHSLGGGTGSGMGTLLISKVREEYPDRIMNTFSVVPSPKVSDTVVE
PYNATLSVHQLVENTDETYCIDNEALYDICFRTLKLAVNMVPFPRNVKEISFVDWCPTGFK
VGINYQPPTVVPGGDLAKVQRAVCMLSNTTAIAEAWARLDHKFDLMYAKRAFVHWYVGEG
MEEGEFSEARED
```
|Name|Description|
|---|----|
|sp| Swiss-Prot|It is a high-quality, manually annotated protein sequence database maintained by the Swiss Institute of Bioinformatics (SIB). |


When a protein sequence is obtained from Swiss-Prot, the sequence identifier in the FASTA format is usually prefixed with "sp|" followed by the accession number of the protein. For example, in the above FASTA format sequence identifier:

"sp" indicates that the protein sequence is from Swiss-Prot, "P69905" is the accession number of the protein, and the rest of the information provides details about the protein's origin, function, and other properties.

Sequences are expected to be represented in the standard IUB/IUPAC amino acid and nucleic acid codes, with these exceptions:
lower-case letters are accepted and are mapped into upper-case;
a single hyphen or dash can be used to represent a gap of indeterminate length;
in amino acid sequences, U and * are acceptable letters (see below).
any numerical digits in the query sequence should either be removed or replaced by appropriate letter codes (e.g., N for unknown nucleic acid residue or X for unknown amino acid residue).

The nucleic acid codes are:
|Code| Nucleic Acid|Code|Nucleic Acid|
|---|----|---|----|
| A --> |adenosine |          M --> |A C (amino)|
|C --> |cytidine    |        S --> |G C (strong)|
|G -->| guanine   |          W --> |A T (weak)|
|T --> |thymidine  |         B -->| G T C|
|U --> |uridine    |         D -->| G A T|
|R -->| G A (purine)  |      H --> |A C T|
|Y -->| T C (pyrimidine) |   V --> |G C A|
|K -->| G T (keto)  |        N --> |A G C T (any)|
                                  -  gap of indeterminate length
                                  
                                  
The accepted amino acid codes are:
|Code| Abbreation|Amino Acid|Code|Abbreation|Amino Acid|
    |---|----|---|----|
    |A |ALA |alanine |                       | P| PRO |proline|
   | B| ASX |aspartate or asparagine |      |  Q |GLN |glutamine|
   | C |CYS |cystine    |                   |  R |ARG |arginine|
   | D |ASP |aspartate  |                   |  S |SER |serine|
  |  E| GLU |glutamate  |                  |   T |THR |threonine|
  |  F |PHE |phenylalanine |                |  U  |   selenocysteine|
  |  G |GLY |glycine   |                   |   V |VAL |valine|
  |  H |HIS |histidine  |                  |   W| TRP |tryptophan|
  |  I |ILE |isoleucine |                  |   Y |TYR| tyrosine|
  |  K| LYS |lysine   |                     |  Z |GLX |glutamate or glutamine|
 |   L |LEU |leucine   |                   |   X  |   any|
 |   M |MET| methionine |                  |   *  |   translation stop|
  |  N |ASN |asparagine   |                |   -   |  gap of indeterminate length|
