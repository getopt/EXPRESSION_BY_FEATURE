### DESCRIPTION

The purpose of **summarize_features.py** is to process the output of
[rpkm_genic_features.bxt5.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/rpkm_genic_features.bxt5.md).
For every gene, **summarize_features.py** summarizes read counts and computes
RPKM in following six types of features:

 - upstream genic regions
 - 3'UTRs
 - CDS
 - introns
 - 5'UTRs
 - downstream genic regions

**summarize_features.py** distinguishes between reads aligned to forward and
reverse genomic strands and it also takes into account gene strand orientation.
In the end, for each gene, the summary read counts and RPKMs are computed separately
for "sense" and "antisense" reads.

**summarize_features.py** unambiguously distributes the reads between
overlapping genes. Evaluation of gene expression by features begins from genome
annotation by
[fix_genic_features.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/fix_genic_features.md).
This annotation is a list of feature coordinates that are strictly
non-overlapping with the respect to an individual gene. However, features of
different genes frequently overlap. **summarize_features.py** distributes read
counts between overlapping genes.

Overlaps between coding genes and non-coding genes are treated by
**summarize_features.py** in a specific manner. Biology of several types of
non-coding RNA genes is special and some of the non-coding RNA transcripts are
frequently extremely abundant in NGS libraries. To avoid artificial elevation
of read counts of coding genes, when overlaps with coding genes are detected
**summarize_features.py** conservatively assigns all "sense" and "antisence"
reads exclusively to the non-coding RNA genes.



### INPUT ARGUMENTS AND FORMAT 

##### -p or --plusTab "path/to/plus strand table" 
##### -m or --minusTab "path/to/minus strand table"

Input to **summarize_features.py** is two tables produced by
[rpkm_genic_features.bxt5.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/rpkm_genic_features.bxt5.md).
Each of the two tables is a list of genomic regions (standard
TAB-delimited genome position format), where each region is annotated with
following info:

1. One of the six types of genic features (e.g. *exon*)
2. Gene name (e.g. *Hsp70Bbb*) 
3. Gene strand orientation (e.g. *+*)
4. Uniquely mappable fraction of the region (e.g. *0.03*)
5. Mappability adjusted RPKM (e.g. *20.183*)
6. Count of reads (e.g. *3*)

When multiple genes overlap with a region, then information in **1**, **2**,
**3** is given for each gene. Exception to these rules are intergenic regions:
they are not associated with any gene specific information and instead are
labelled simply `intergenic`.

Two input tables have to be supplied to **summarize_features.py**: one is for
reads aligned to forward (or plus) genomic strand, and another table is for
reads aligned to the reverse (or minus) strand. Currently it is required for
the two tables to be of the same length and to list genomic regions in the same
order.

##### -n or --mappedReads "number" *or* "path/to/bowtie1 log file"

It is required to supply the number of mapped reads in order to compute RPKM.
Reads can be supplied on the command line by giving the actual number (e.g.
'10000000') or by supplying a path to bowtie1 log file (i.e. file with standard
error output of bowtie1). Such files contain the number of maped reads and
**summarize_features.py** parses it out.

##### --tableType 'gene' or 'partition'

This argument when set to 'gene' instructs to output the table with one gene
per row and RPKMs and read counts per features of that gene. This is the
default mode. Setting the argument to 'partition' instructs to output the table
per one partition type where counts from such partitions in from all genes are
summarized by one number. The 'partition' mode is useful to get very zoomed-out
view of the library, for example, to see whether expression of `utr3p` features
is higher than `utr5p` features on a global scale. The 'partition' mode is
currently under development. 


### PROCEDURE IN DEFAULT MODE (--tableType 'gene')

**summarize_features.py** processes one line at a time from *plus strand table*
and *minus strand table* input files. Each line in the input correspons to a
genomic region that overlaps with features of at least one gene: either
`utr5p`, `exon` and `utr3p` features or `upstream` and `downstream` features
(purely `intergenic` regions are skipped in --tableType 'gene' mode).

To correctly distribute reads mapped to the region, the program first parses
out names of genes whose features overlap in the region. Then reads mapped to
the region are distributed evenly between overlapping genes, taking
strandedness of genes and of read alignments into account. Whenevery possible,
reads are counted as "sense". Reads are counted as "antisense" only when there
is no genes in sense orientation overlapping the region. Note, that regions
that overlap features of non-coding and/or repetitive genes are treated
specially. These principles are illustrated by the below examples:


