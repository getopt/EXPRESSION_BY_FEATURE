### DESCRIPTION [bin/fix_genic_features.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/bin/fix_genic_features.py)

### PROCEDURE

`fix_genic_features.py` makes in total three passes through annotation to
arrive at the final genome wide annotation.

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
  recieve final `upstream` and `downstream` annotation can be shorter than 1kb:
  if we detect that the upstream region of geneA overlaps with `utr5p`,
  `utr3p`, `exon` or `intron` of geneB, then only the non-overlaping fraction
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

##### -t or --gtf "path/to/transcript annotation in GTF fromat"

Input transcript annotation is required to be in GTF format
(http://genome.ucsc.edu/FAQ/FAQformat.html#format4). GTF RefSeq transcript
annotation can be obtained via UCSC Table Browser
(http://genome.ucsc.edu/cgi-bin/hgTables?command=start). Unfortunately the
default GTF table does not include `gene_name` attributes. Fortunately,
however, the relationship between RefSeq `transcript_id` and `gene_name` field
is preserved when UCSC Table browser is instructed to output `all fields`.
Therefore it easy to supplement the GTF table with `gene_name` attributes via
Perl snippet
[add_gene_names_to_gtf.pl](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/bin/add_gene_names_to_gtf.pl). 

Example input GTF file with added `gene_name` attributes:

```
...
chr4    dm3_refGene exon        1048489 1048508 0.000000    +   .   gene_id "NR_073624"; transcript_id "NR_073624"; gene_name "CR44031" 
chr4    dm3_refGene exon        1049336 1049985 0.000000    +   .   gene_id "NR_073624"; transcript_id "NR_073624"; gene_name "CR44031" 
chr4    dm3_refGene stop_codon  892119  892121  0.000000    -   .   gene_id "NM_143692"; transcript_id "NM_143692"; gene_name "unc-13" 
...
```

Note: only `exon`, `stop_codon` and `start_codon` features are important. Other
features (such as `CDS`) can be present in the GTF input file, but they are
ignored by `fix_genic_features.py`.

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
with chromsome names and lengths can be generated with
[get_chr_lengths.pl](https://github.com/getopt/FASTA_TOOLS/blob/master/get_chr_lengths.pl)


##### -a or --assembly "name of genome assembly"

Possible values here are *dm3*, *mm10* and *hg19*. 

We would like to recognize and label non-coding genes (e.g. miRNAs, snoRNAs,
CR-genes, etc.) at an early stage. To do this, during processing GTF
annotation, `fix_genic_features.py` identifies gene names that belong to
categories of non-coding RNA genes (and also histone genes in case of *dm3*)
and joins `gene_name` and `transcript_id` attributes via a `#` symbol for
several selected types of genes. Introduction of `#` symbol is convenient
during subsequent steps of analysis as a lable of genomic partitions that
overlap with non-coding genes. Identification of non-coding RNA genes is done
via pattern matching, which is dependent on the type of the organizm in
question (since gene naming schema differs between mouse, human and fly).
Currently following flags are used to recognized non-coding genes:

```
dm3  : 'mir-' '^CR'  ':' 
mm10 : 'Mir'  'Snor' 'Rik$'
hg19 : 'MIR'  'SNOR'
```

Note:
- Matching by ':' in *dm3* also matches histone genes.
- Joining `genen_name` `transcrip_id` attributes is very important for
  gene-name based annotation for some categories of genes (e.g. miRNA genes).
  Some genes do not always have unique `gene_name` attributes despite location
  in distinct genomic locations, and for them attaching `transcript_id` to
  `gene_name` is a way to create a unique label of each genomic site.
- During the first pass of `fix_genic_features.py`, every genic feature on any
  non-coding gene (i.e. `exon` and `intron`) is expanded by 25 nt from both
  ends. This extension of features of non-coding genes is done since during
  subsequent counting stage (counting is performed by
  [rpkm_genic_features.bxt5.py](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/rpkm_genic_features.bxt5.md))
  we really do not want to count any of the reads overlapping with non-coding
  genes as if they belong to a coding gene. Therefore we artificially expand
  coordinates of features of non-coding genes by 25 nt. Note, that even though
  both `exon` and `intron` features are originally expanded, only `exon`
  features remain expanded in the final output, since `exon` trumps `intron`
  during the third pass of `fix_genic_features.py` when ambiguities of
  overlapping features are resolved (see `PROCEDURE` section above for more
  details). 
