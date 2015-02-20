### DESCRIPTION [bin/fix_genic_features.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/bin/fix_genic_features.py)

### PROCEDURE

`fix_genic_features.py` makes in total three passes through transcript
annotation to arrive at the final genome wide annotation.

In the first pass, for every distinct `gene_name` attribute in GTF file we
store all of the corresponding features from the GTF input file (i.e. all
exons, and also all start and stop codons in case of coding genes).

In the second pass, we define the six types of genomic features. Features can
be overlapping at first.

 - upstream genic regions
 - 3'UTRs
 - CDS
 - introns
 - 5'UTRs
 - downstream genic regions

In the third and final pass, we go through the list of features and resolve
ambiguities of overlapping features using following principles:

- Any exon outside of outermost start codon is labelled `utr5p`, outside of the
  outermost stop codon is labelled `utr3p` and middle exons are labelled simply
  `exon`.

- If exon overlaps with an intron, we label such region `exon`.

- Upstream and downstream regions of any gene are initially defined in the
  second pass as 1 kb up- and downstream. During the third pass regions that
  receive final `upstream` and `downstream` annotation can be shorter than 1kb:
  if we detect that the upstream region of geneA overlaps with `utr5p`,
  `utr3p`, `exon` or `intron` of geneB, then only the non-overlapping fraction
  of the upstream regions of geneA is going to be labelled `upstream` in the
  final annotation. In the extreme cases the `upstream` and `downstream`
  recieve length of 0 nt. For example, this can happen when geneA is located
  within a long intron of geneB and the intial 1 kb upstream and downstream
  regions of geneA happen to be entirely contained within the long intron of
  geneB.


##### Definition of exonic and intronic regions

Exons are defined in the first pass from the information in the input GTF file.
Regions between consecutive exons are all labelled as `intron` in the second
pass. Some exons will be re-labelled into `utr5p` and `utr3p` in the second
pass. Note regions of introns (or sometimes entire introns) lose `intron` label
in the third pass if an overlap with exonic features is detected.

##### Definition of UTRs

UTR features of coding genes are defined in the second pass after we store
coordinates of all exons, start and stop codons. In the second pass we select a
single start and a single stop codon position per gene. Then, these two
selected positions are used as stable reference points for assigning UTR
annotation to `exon` features of all transcripts of that gene.  When gene
encodes multiple transcript isoforms that differ in positions of start and stop
codons, we choose the positions of the *outermost* codons as the reference
points. Exons, or parts of exons, that are outside of the reference points
receive `utr5p` and `utr3p` labels, while middle exons remain with `exon`
label. In case of non-coding genes, i.e. genes without `start_codon` and
`stop_codon` entries in the input GTF file, all exons remain with the `exon`
label.

##### Definition of upstream and downstream regions

Upstream and downstream regions are defined as 1kb regions outside of outermost
exons in the second pass of the program. They will be contracted in the third
pass if overlaps with genic features is detected.

### INPUT ARGUMENTS AND FORMAT 

Input data is passed to `fix_genic_features.py` via three arguments.  Command
line flags for passing these arguments to the program, and description of input
format is given below.

##### -g or --genomeFile "path/to/file with chromosome names and lengths"

To correctly annotate regions outside of genes (i.e. regions that recieve
`intergenic` label) it is required to know chromosome lengths. This file 
consists of lines with two tab-delimited fields. First field is chromsome
name, the second field is chromsome length. 

Example of input genome file for *Drosophila melanogaster* dm3 genome assembly:

```
chr2L       23011544
chr2LHet    368872
chr2R       21146708
chr2RHet    3288761
chr3L       24543557
chr3LHet    2555491
chr3R       27905053
chr3RHet    2517507
chr4        1351857
chrM        19517
chrU        10049037
chrUextra   29004656
chrX        22422827
chrXHet     204112
chrYHet     347038
```
Note: given that we have the genome sequence file in FASTA format, the table
with chromosome names and lengths can be generated with
[get_chr_lengths.pl](https://github.com/getopt/FASTA_TOOLS/blob/master/get_chr_lengths.pl)

##### -t or --gtf "path/to/transcript annotation in GTF fromat"

Input transcript annotation is required to be in GTF format
(http://genome.ucsc.edu/FAQ/FAQformat.html#format4). GTF RefSeq transcript
annotation can be obtained via UCSC Table Browser
(http://genome.ucsc.edu/cgi-bin/hgTables?command=start). 

In the GTF file only `exon`, `stop_codon` and `start_codon` features are
important. Other features (such as `CDS`) can be present in the GTF input file,
but they are ignored by `fix_genic_features.py`. Things like intron are
inferred based `exon` position.

For it's function of per-gene annotation of genomes, `fix_genic_features.py`
requires presence of the `gene_name` attribute in the input GTF file, since
annotation is done on a per-gene basis. Unfortunately the default GTF table
does not include `gene_name` attributes. Fortunately, however, the
relationship between RefSeq `transcript_id` and `gene_name` field is preserved
when UCSC Table browser is instructed to output `all fields`.  Therefore it
easy to supplement the GTF table with `gene_name` attributes via Perl snippet
[add_gene_names_to_gtf.pl](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/bin/add_gene_names_to_gtf.pl). 

In later stages of the analysis, performed by
[summarize_features.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/bin/summarize_features.py)
we distribute reads between overlapping genes. Overlaps of two protein coding
genes are treated differently from overlaps of protein-coding and non-coding
genes (this is because some highly abundant non-coding RNAs are artifacts of
cloning see [documentation for
summarize_features](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/summarize_features.md)).
Therefore we manually curate GTF files from UCSC Table browser and add
`gene_type` attribute to each feature (e.g. `protein_coding`,
`AR_miRBase_hairpin`, and etc.). The curation of the UCSC RefSeq GTF files is
described
[here](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/curating_RefSeq_GTF.md).

Example input GTF file with added `gene_name` and `gene_type` attributes:

```
...
chr1    mm10_refGene    exon    3421702 3421901 0.000000    -   .   gene_id "NM_001011874"; transcript_id "NM_001011874"; gene_name "Xkr4"  ; gene_type "protein_coding";
chr1    00RM00000000    exon    3428386 3428449 .   +   .   gene_id "NR_RM:SRP"; transcript_id "NR_RM:SRP"; gene_name "NR_RM:7SLRNA_dup1"; gene_type "AR_SRP";
chr1    00RM00000000    exon    3529769 3529806 .   +   .   gene_id "NR_RM:tRNA"; transcript_id "NR_RM:tRNA"; gene_name "NR_RM:tRNA-Gly-GGA_dup1"; gene_type "AR_tRNA";
chr1    mm10_refGene    CDS 3670552 3671348 0.000000    -   0   gene_id "NM_001011874"; transcript_id "NM_001011874"; gene_name "Xkr4"  ; gene_type "protein_coding";
...
```

The `gene_type` attribute is concatenated to `gene_name` attribute via `#`
character in the output file. 


##### How we make sure that gene names are unique

Since `fix_genic_features.py` uses `gene_name` attribute as a key in a
dictionary to store exon information, it is a problem if several genes share
the same name. Luckily, the UCSC RefSeq GTF file contains `transcript_id`
attribute that is strictly unique per-chromosome: when copies of the same
sequence present as different genes on one chromosome in the GTF file they are
distinguished by `_dupN` suffix. Therefore, we reasoned that the easiest way
to make sure that `gene_name` attributes are unique, `fix_genic_features.py`
simply moves `_dupN` suffix from `transcript_id` onto `gene_name` whenever the
suffix is present.

When we manually curate the GTF file and add some non-coding transcripts to it
(description of curating procedure is
[here](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/curating_RefSeq_GTF.md)),
we make sure that `gene_name` is unique to start with by adding `_dupN` suffx
to the `gene_name` straight away.


##### Conservative extension of gene boundaries belonging to certain `gene_type` categories

During the first pass of `fix_genic_features.py`, every genic feature on any
non-coding gene belonging to selected `gene_type` categories  is expanded by 25
nt from both ends. Currently this procedure is performed for `AR_miRBase_hairpin`,
`AR_snRNA`, `AR_snoRNA`.

This extension of features of non-coding genes is done since during subsequent
counting stage (counting is performed by
[rpkm_genic_features.bxt5.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/rpkm_genic_features.bxt5.md))
we really do not want to count any of the reads overlapping with non-coding
genes as if they belong to a coding gene. Therefore we artificially expand
coordinates of features of non-coding genes by 25 nt. Note, that even though
both `exon` and `intron` may originally expanded, only `exon` features remain
expanded in the final output, since `exon` trumps `intron` during the third
pass of `fix_genic_features.py` when ambiguities of overlapping features are
resolved (see `PROCEDURE` section above for more details). 
