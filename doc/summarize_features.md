### DESCRIPTION

The purpose of `summarize_features.py` is to process the output of
[rpkm_genic_features.bxt5.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/rpkm_genic_features.bxt5.md).
For every gene, `summarize_features.py` summarizes read counts and computes
RPKM in the following six types of features:

 - upstream genic regions
 - 3'UTRs
 - CDS
 - introns
 - 5'UTRs
 - downstream genic regions

`summarize_features.py` distinguishes between reads aligned to forward and
reverse genomic strands and it also takes into account gene strand orientation.
In the end, for each gene, the summary read counts and mappability adujusted
RPKMs are computed separately for "sense" and "antisense" reads.

`summarize_features.py` unambiguously distributes the reads between
overlapping genes. Evaluation of gene expression by features begins from genome
annotation by
[fix_genic_features.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/fix_genic_features.md).
This annotation is a list of feature coordinates that are strictly
non-overlapping with the respect to an individual gene. However, features of
different genes frequently overlap. `summarize_features.py` distributes read
counts between overlapping genes.

Overlaps between coding genes and non-coding genes are treated by
`summarize_features.py` in a specific manner. Biology of several types of
non-coding RNA genes is special and some of the non-coding RNA transcripts are
frequently extremely abundant in NGS libraries. To avoid artificial elevation
of read counts of coding genes, when overlaps with coding genes are detected
`summarize_features.py` conservatively assigns all "sense" and "antisence"
reads exclusively to the non-coding RNA genes.


### PROCEDURE IN DEFAULT MODE (--tableType 'gene')

`summarize_features.py` processes two input files *plus strand table* and
*minus strand table* one line at a time. Each line in the input describes a
genomic region that overlap features of one or more genes. These features are
either `utr5p`, `exon` and `utr3p` or `upstream` and `downstream`. Note that
purely `intergenic` regions are skipped in --tableType 'gene' mode.

In order to distribute reads mapped to a region described by a line in the
input file, `summarize_features.py` first parses out names of one or more genes
whose features overlap the region. Then reads mapped to the region are
distributed evenly between the genes, taking strandedness of genes and
of read alignments into account. If the region containes a feature of a single
gene, then reads will be distributed as "sense" or "antisense" with the respect
to the gene strandedness. When multiple genes overlap the region, reads are
distributed to be in "sense" orientation whever possible. Principles of
`summarize_features.py` are illustrated by the following examples:

1. A region overlaps `intron` feature of gene_A that is on (+)-strand and
   `exon` feature of gene_B that is also on (+)-strand. There are 10 reads aligned
   to the (+)-strand, and 2 reads aligned to the (-)-strand. The program
   distributes counts of the aligned reads as:

    - 5 "sense" reads to `intron` of gene_A
    - 5 "sense" reads to `exon` of gene_B
    - 1 "antisense" read to the `intron` of gene_A
    - 1 "antisense" read to the `exon` of gene_B

2. If in the example **1** the strandedness of gene_B is chaged to (-)-strand,
   then the outcome becomes:

    - 10 "sense" reads to `intron` of gene_A
    - 2  "sense" reads to `exon` of gene_B  

3. Now the region overlaps `intron` feature of gene_A and `exon` feature of
   gene_B, both on (+)-strand, and also `exon` of gene_C on (-)-strand. As in
   **1** and **2**, there are 12 reads mapped to the region (10 to the (+)-strand,
   and 2 to (-)-strand). In this case:
    
    - 5 "sense" reads to `intron` of gene_A
    - 5 "sense" reads to `exon` of gene_B
    - 2 "sense" reads to `exon` of gene_C

4. Special case is when a region overlaps features of one or more non-coding
   repetitive genes. In such cases, all read counts are distributed to the
   non-coding gene or genes, while 0 reads is attributed to the coding genes. For
   example, if gene_A, gene_B and gene_C from case *4* are coding, but the
   region also overlaps `exon` of non_coding_gene_A on (+)-strand and `exon` of
   non_coding_gene_B on (-)-strand, then:
    
    - 10 "sense" reads to `exon` of non_coding_gene_A
    - 2 "sense" reads to `exon` of non_coding_gene_B
    - 0 reads to features of gene_A, gene_B and gene_C


### INPUT ARGUMENTS AND FORMAT 

##### -p or --plusTab "path/to/plus strand table" 
##### -m or --minusTab "path/to/minus strand table"

Input to `summarize_features.py` is two tables produced by
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

When multiple genes overlap with a region, then the input files contain **1**,
**2**, **3** information for each gene. Exception to these rules are intergenic
regions: they are not associated with any gene specific information and instead
are labelled simply `intergenic`.

Two input tables have to be supplied to `summarize_features.py`: one is for
reads aligned to forward (or plus) genomic strand (*plus strand table*), and
another table is for reads aligned to the reverse (or minus) strand (*minus
strand table*). Currently it is required for the two tables to be of the same
length and to list genomic regions in the same order.

##### -n or --mappedReads "number" *or* "path/to/bowtie1 log file"

It is required to supply the number of mapped reads in order to compute
mappability adjusted RPKM.  Reads can be supplied on the command line by giving
the actual number (e.g.  '10000000') or by supplying a path to bowtie1 log file
(i.e. file with standard error output of bowtie1). Such files contain the
number of maped reads and `summarize_features.py` parses it out.

##### --tableType 'gene' or 'partition'

This argument when set to 'gene' instructs to output the table with one gene
per row and mappability adjusted RPKMs and read counts per features of that
gene. Setting the argument to 'partition' instructs to output the table per one
partition type where counts from such partitions in from all genes are
summarized by one number. The 'partition' mode is useful to get very zoomed-out
view of the library, for example, to see whether expression of `utr3p` features
is higher than `utr5p` features on a global scale. Currently 'gene' mode is
default, while the 'partition' mode is under development. 

