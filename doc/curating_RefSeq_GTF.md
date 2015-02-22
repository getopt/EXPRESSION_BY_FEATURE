### DESCRIPTION OF HOW REFSEQ GTF FILES ARE CURATED**


*work in progress*

Following attributes are added to each feature of genes in the GTF file:

```
AR_RNA7K
AR_SRP
AR_miRBase_hairpin
AR_rRNA
AR_scRNA
AR_snRNA
AR_snoRNA
AR_tRNA
lncRNA_encoding_wkncRNA
ncRNA_other
protein_coding
```

*work in progress*

RefSeq contains many non-coding RNA annotations (feature names starting with NR_). However, it is difficult to infer the type of RNA by its name, as for examples some long non-coding transcripts are named 'mir' for microRNAs. This is undesirable, as in donwstream analysis microRNAs, tRNAs, etc. overlapping with protein-coding genes get all reads, while the corresponding overlapping protein-coding gene regions get none. Thus, long pseudo-microRNAs can "absob" reads mapping to large regions.
To overcome this, we gather non-coding RNA annotations from different databases, including

```
RFAM (bigBed format, 0-based coordinates); abbreviated RF; we take annotations of all well known RNAs excluding microRNAs
Ensembl (gtf format, 1-based coordinates); abbreviated EN; we take annotations of all well known RNAs excluding microRNAs
RepeatMasker (table, 1-based coordinates); abbreviated RM; we take annotations of all well known RNAs excluding microRNAs
miRBase (gff3 format, 1 based coordinates); abbreviated MB; for microRNA annotations only
flyBase (gff format, 1 based coordinates); abbreviated FB; we take annotations of all well known RNAs excluding microRNAs
```
Well-known non-coding RNAs of the classes outline above are filtered from each database by their name using the following patterns which were established by manual inspection of individual files:

```
RNA7K='\b(7SK|Rn7sk)'
tRNA='tRNA'
snRNA='\b(snRNA|U[12]\b|U[4567]a?|U1[12])'
snoRNA='\b(scaRNA|snoRNA|ACEA\_U3|ACA64|sno(R|U|Z)*|SNOR(A|D|RD)*|U[38])'
rRNA='rRNA'
scRNA='scRNA'
SRP='(SRP|srpRNA)'
```
and a temporary file with the following naming is created: "wkncRNA.tRNA.tmp";

To remove redundancies, we then use the bedtools command 'merge' and create a file where for each overlap we have the end coordinates, and all names. Example

```
RFAM table:
chr1	RF	exon	7265804	7265932	.	.	.	SNORA17 (note there is no strand info)

Ensembl table:
chr1	EN	exon	7265803	7265932	.	+	.	gene_id&"ENSMUSG00000077244";&gene_version&"1";&gene_name&"Gm23274";&gene_source&"EN";&gene_biotype&"snoRNA";

Resulting table:
chr1	EN;RF	exon	13774333	13774464	.	-;.	.	gene_id&"ENSMUSG00000077368";&gene_version&"1";&gene_name&"Gm26273";&gene_source&"EN";&gene_biotype&"snoRNA";|SNORA17
```

These names and sources of coordinates are then parsed into the file wkncRNA.snoRNA.tmp.merged.parsed:

```
chr1	RF00EN000000	exon	7265802	7265932	.	+	.	gene_id "NR_RFEN:snoRNA"; transcript_id "NR_RFEN:snoRNA"; gene_name "NR_EN:Gm23274,RF:SNORA17_dup1"; gene_type "AR_snoRNA";
```
The second column shows all databases that were used to annotate this region as a gene\_type snoRNA; gene\_id and transcript\_id fields are created based on the sources and type (RFEN:snoRNA means snoRNA supported by RF and EN annotations); the gene\_name is a comma separated list of names in different databases. To avoid duplicated names, the first instance of each name has a 'dup1' tag, and ech next instance of the same name has 'dup2', 'dup3', etc. tag. 
For strandedness, if all sources consistently show minus or plus strand, this is peserved. If some sources show plus (or minus) strand, and others have no strand info, then we assume the strand info of the source that has it. If there is conflicting strandedness, we print two lines, one with positive and one with negative strand.
```
+;+ becomes +
+;. becomes +
-;. becomes -
-;- becomes -
+;- creates two lines, one with + and one with -
+;-;. same as above
```
we then put together all different wkncRNA files into one big wkncRNAs.gtf named after the species genome, i.e. mm10\_wkncRNAs.gtf, dm3\_wkncRNAs.gtf and hg19\_wkncRNAs.gtf.







