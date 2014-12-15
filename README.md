#####EXPRESSION_BY_FEATURE REPOSITORY DESCRIPTION


### RATIONALE 

A single gene frequently encodes multiple overlapping transcripts and genes
themselves occasionally overlap. Therefore classification of genomic regions
into exonic, intronic and other related categories can be ambiguous.

**fix_genic_features.py** examines transcript annotation and resolves ambiguity
of genomic annotation on per-gene basis. The program uses a hierarchy principle
to resolve ambiguities, that for most protein coding genes allows to eventually
attribute following six types of non-overlapping genomic features:

 - upstream genic regions
 - 3'UTRs
 - CDS
 - introns
 - 5'UTRs
 - downstream genic regions

The output of **fix_genic_features.py** is a whole-genome annotation. This
annotation allows to quickly characterize every unique alignment of an NGS
library with **rpkm_genic_features.bxt5.py** program, and then to quantify expression
of every gene and of every feature-type by **summarize_features.py**. 


### DOCUMENTATION 

Detailed documentation of programs of this repository is in doc/ directory:

[fix_genic_features.py documenation](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/fix_genic_features.md)

[rpkm_genic_features.bxt5.py documentation](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/rpkm_genic_features.bxt5.md)

[summarize_features.py documentation](https://github.com/getopt/EXPRESSION_BY_FEATURE/blob/master/doc/summarize_features.md)



