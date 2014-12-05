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

Evaluation of gene expression by features begins from genome annotation by
[fix_genic_features.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/fix_genic_features.md).
This annotation is a list of feature coordinates that are strictly
non-overlapping with the respect to an individual gene. However, features of
different genes frequently overlap. In such cases **summarize_features.py**
unambiguously distributes the reads between overlapping genes. 

Overlaps between coding genes and non-coding genes are treated by
**summarize_features.py** in a specific manner. Biology of several types of
non-coding RNA genes is special and some of the non-coding RNA transcripts are
frequently extremely abundant in NGS. To avoid artificially elevated read
counts of coding genes that overlap with non-coding RNA genes,
**summarize_features.py** conservatively assigns reads from such overlapping
regions to non-coding RNA gene only and never to the coding gene (irrespective
of "sense" or "antisense" orientation).


### INPUT ARGUMENTS AND FORMAT 

##### -p or --plusTab "path/to/plus strand table" *and* -m or --minusTab "path/to/minus strand table"

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

Two input tables have to supplied to **summarize_features.py**: one is for
reads aligned to forward (or plus) genomic strand, and another table is for
reads aligned to the reverse (or minus) strand. It is required for the two
tables be of the same length and to list genomic regions in the same order.

##### -n or --mappedReads "number" *or* "path/to/bowtie1 log file"

It is required to supply the number of mapped reads in order to compute RPKM.
Reads can be supplied on the command line by giving the actual number (e.g.
'10000000') or by supplying a path to bowtie1 log file (i.e. file with standard
error output of bowtie1). These log files are of standard format, and
**summarize_features.py** parse out the number of mapped reads from it.

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

**summarize_features.py** 
