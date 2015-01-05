#!/usr/bin/env python

''' Usage:

-p or --plusTab          <tab>  path/to/table_for_(+)_genomic_strand
-m or --minusTab         <tab>  path/to/table_for_(+)_genomic_strand

-n or --mappedReads      <int>  total number of mapped reads i.e. number that
                                we want to use for RPKM calculation

                                HACK FOR PIPELINE ENVIRONMENT: instead of the
                                number, a path to bowtie log file can be
                                supplied and the number will be parsed out from
                                it this is handy in the pipeline, when number
                                of mapped reads is unknown at the time of
                                creating job submission

-t or --type  <partition|gene>  determine the type of table that will be
                                produced. Set to 'gene' to make table for
                                individual genes.  Set to 'partition' to make
                                table for large categories, like 'exon',
                                'intron'. DEFAULT: 'gene'


INPUT: two tables (one for reads mapped to the forward [plus] genomic strand,
       and one for reads mapped to the reverse [minus] genomi strand, in
       '*.per_geneSS_[plus|minus].tab' format, looks like:

...
chr2L 102086 102379 downstream~Zir{+};upstream~CG11377{+}  1.0  0     0.0
chr2L 102379 102458 utr5p~CG11377{+}                       1.0  0     0.0
chr2L 102458 102906 exon~CG11377{+}                        1.0  49.5  86.0
chr2L 102906 103005 intron~CG11377{+}                      1.0  0     0.0
chr2L 103005 103434 exon~CG11377{+}                        1.0  13.8  23.0
chr2L 103434 103515 intron~CG11377{+}                      1.0  0     0.0
chr2L 103515 103936 exon~CG11377{+}                        1.0  14.0  23.0
chr2L 103936 103961 utr3p~CG11377{+}                       1.0  51.5  5.0
chr2L 103961 104142 utr3p~CG11377{+};utr3p~Nhe1{-}         1.0  0     0.0
...

Where columns are:
    1. Chromosome name
    2. Start (0-based)
    3. End (1-based)
    4. Annotation of the region 
    5. Fraction of the region mappable by reads with unique alignments
    6. RPKM adjusted for mappability
    7. Raw read count



TODO:

   - Add summarize_per_partition()

   - Modify summarize_per_gene() to work on one input file at a time (e.g.
     first on --plusTab <table> and then on --minusTab <table>, rather then on
     two files at the same times). We don't actually have to check that they
     have the same number of lines and that they are in the same order.
     Removing this constrain is going to allow to easily work with
     not-strand-specific regular RNA-seq data.

'''

from __future__ import print_function
from __future__ import division
import pysam
import getopt
import sys
import re
import collections

# hardcoded, because only one type of tables is parsed by this script
CHROMCOL = 0
STARTCOL = 1
ENDCOL   = 2
NAMECOL  = 3
FRACCOL  = 4
RPKMCOL  = 5
COUNTCOL = 6
PARTITIONS_FOR_GENES = \
['upstream', 'utr5p', 'exon', 'intron', 'utr3p', 'downstream']


def getOptions( argv ):
    plusTabfile   = ''
    minusTabfile  = ''
    mappedReads   = ''
    tableType     = 'gene'
    try:
        opts, args = getopt.getopt(argv[1:], "hp:m:n:t:", \
        [ "help", "plusTab=", "minusTab=", "mappedReads=", "type=" ])
    except getopt.GetoptError:
        print ('Usage: %s -p <plusTab> -m <minusTab> -n <mappedReads> \
                                      -t <type>'  % argv[0], file = sys.stderr) 
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print ('Usage: %s -p <plusTab> -m <minusTab> -n <mappedReads> \
                                       -t <type>' % argv[0], file = sys.stderr)
            sys.exit()
        elif opt in ("-p", "--plusTab"):
            plusTabfile = arg
        elif opt in ("-m", "--minusTab"):
            minusTabfile = arg
        elif opt in ("-n", "--mappedReads"):
            mappedReads = arg
        elif opt in ("-t", "--type"):
            tableType = arg
        else:
            assert False, "unhandled option"
    return plusTabfile, minusTabfile, mappedReads, tableType


def parse_mappedReads( something ):
    '''
    This is a hack that allows to supply with -m (i.e. --mappedReads) either
    the number of mapped reads or a path to bowtie error log file. The log file
    name has to include "bowtie" in its name. If the path to the log file is
    supplied, then mapped reads are being parsed out from the content of the
    file. Otherwise, the int(argument) is returned.
    '''
    # check whether something supplied with --mappedReads 
    # looks like bowtie error log file
    if "bowtie" in something:                                                       
        bowtieLogHandle = open(something, 'r')
        for line in bowtieLogHandle:
            # found the line that contains info about mapped reads
            if "reads with at least one reported alignment:" in line:               
                m = re.search(':\s(\d+)\s\(', line)
                # NB: in python group(0) returns entire match; 
                # group 1 corresponds to perl $1
                return int(m.group(1))                                              
    else:                                                                         
        # else means that something supplied with --mappedReads is 
        # presumably a number of mapped reads
        return int(something)


def parse_line ( line ):
    fields     =  line.strip().split('\t')
    chrom      =  fields[CHROMCOL]
    start      =  int(fields[STARTCOL])
    end        =  int(fields[ENDCOL])
    name       =  fields[NAMECOL]
    frac       =  float(fields[FRACCOL])
    rpkm       =  float(fields[RPKMCOL])
    count      =  float(fields[COUNTCOL])
    length     =  end - start
    mapLength  =  length * frac
    return chrom, start, end, name, rpkm, count, length, mapLength


def initialize_geneID( summary, geneID ):
    '''
    For every newly encountered gene, initialize an entry in summary{}
    '''
    summary[geneID] = {}

    summary[geneID]['chrom']  = ''
    summary[geneID]['start']  = ''
    summary[geneID]['end']    = ''
    summary[geneID]['strand'] = ''

    summary[geneID]['sense'] = {}
    summary[geneID]['sense']['upstream']                    = {}
    summary[geneID]['sense']['upstream']["MAPLENGTH"]       = 0
    summary[geneID]['sense']['upstream']["COUNT"]           = 0
    summary[geneID]['sense']['utr5p']                       = {}
    summary[geneID]['sense']['utr5p']["MAPLENGTH"]          = 0
    summary[geneID]['sense']['utr5p']["COUNT"]              = 0
    summary[geneID]['sense']['exon']                        = {}
    summary[geneID]['sense']['exon']["MAPLENGTH"]           = 0
    summary[geneID]['sense']['exon']["COUNT"]               = 0
    summary[geneID]['sense']['intron']                      = {}
    summary[geneID]['sense']['intron']["MAPLENGTH"]         = 0
    summary[geneID]['sense']['intron']["COUNT"]             = 0
    summary[geneID]['sense']['utr3p']                       = {}
    summary[geneID]['sense']['utr3p']["MAPLENGTH"]          = 0
    summary[geneID]['sense']['utr3p']["COUNT"]              = 0
    summary[geneID]['sense']['downstream']                  = {}
    summary[geneID]['sense']['downstream']["MAPLENGTH"]     = 0
    summary[geneID]['sense']['downstream']["COUNT"]         = 0
    summary[geneID]['antisense'] = {}
    summary[geneID]['antisense']['upstream']                = {}
    summary[geneID]['antisense']['upstream']["MAPLENGTH"]   = 0
    summary[geneID]['antisense']['upstream']["COUNT"]       = 0
    summary[geneID]['antisense']['utr5p']                   = {}
    summary[geneID]['antisense']['utr5p']["MAPLENGTH"]      = 0
    summary[geneID]['antisense']['utr5p']["COUNT"]          = 0
    summary[geneID]['antisense']['exon']                    = {}
    summary[geneID]['antisense']['exon']["MAPLENGTH"]       = 0
    summary[geneID]['antisense']['exon']["COUNT"]           = 0
    summary[geneID]['antisense']['intron']                  = {}
    summary[geneID]['antisense']['intron']["MAPLENGTH"]     = 0
    summary[geneID]['antisense']['intron']["COUNT"]         = 0
    summary[geneID]['antisense']['utr3p']                   = {}
    summary[geneID]['antisense']['utr3p']["MAPLENGTH"]      = 0
    summary[geneID]['antisense']['utr3p']["COUNT"]          = 0
    summary[geneID]['antisense']['downstream']              = {}
    summary[geneID]['antisense']['downstream']["MAPLENGTH"] = 0
    summary[geneID]['antisense']['downstream']["COUNT"]     = 0

    summary[geneID]['annotation'] = ''
    summary[geneID]['annotation'] = interprete_geneID( geneID )

    return summary


def interprete_geneID( geneID ):
    # LABEL genes with "#" in the name field: i.e. 'mir-' and 'snoRNA',
    # because '#' is added for special non-coding genes during
    # annoation
    # see
    # "EXPRESSION_BY_FEATURE/fix_genic_features.py"
    annotation = 'NA'
    # drosophila stuff
    if re.search('CR', geneID) and re.search('#', geneID):
        annotation = 'noncoding:CR'
    if re.search('snRNA:', geneID) and re.search('#', geneID):
        annotation = 'noncoding:snRNA'
    elif re.search('snoRNA:', geneID) and re.search('#', geneID):
        annotation = 'noncoding:snoRNA'
    elif re.search('mir-', geneID) and re.search('#', geneID):
        annotation = 'noncoding:miRNA'
    elif re.search('snmRNA:', geneID) and re.search('#', geneID):
        annotation = 'noncoding:snmRNA'
    elif re.search('scaRNA:', geneID) and re.search('#', geneID):
        annotation = 'noncoding:scaRNA'
    elif re.search('rRNA:', geneID) and re.search('#', geneID):
        annotation = 'noncoding:rRNA'
    elif re.search('His', geneID) and re.search('#', geneID):
        annotation = 'histone'
    # mouse stuff
    elif re.search('Rik', geneID) and re.search('#', geneID):
        annotation = 'RIKEN'
    elif re.search('^Mir', geneID) and re.search('#', geneID):
        annotation = 'noncoding:miRNA'
    elif re.search('^Sn', geneID) and re.search('#', geneID):
        annotation = 'noncoding:sn'
    # human stuff
    elif re.search('^MIR', geneID) and re.search('#', geneID):
        annotation = 'noncoding:miRNA'
    elif re.search('^SNOR', geneID) and re.search('#', geneID):
        annotation = 'noncoding:sn'
    # other
    elif re.search('#', geneID):
        annotation = 'noncoding:other'
    else:
        annotation = 'coding'
    return annotation


def gene_entry( summary, geneID, chrom, start, end, geneStrand, partType, mapLength, pCount, mCount, nUnits_plus_str, nUnits_minus_str ):
    '''
    RATIONALE:
    We want to examine a feature of a gene here, and give the gene a fair share
    of "sense" and "antisense" reads for that feature and put it into summary{}
    
    INPUT:
    Therefore input is summary{} and the gene name, togehter with all
    the information about the region and the feature:

     Gene name to which this feature is attributed:
        'geneID'
     
     Strandedness of the gene:
        'geneStrand'
     
     Type of feature in the region:
        'partType'

     Count of reads aligned to + and - strands in the region:
        'pCount' 'mCount'

     Count of genes on + and - strands in the region:
        'nUnits_plus_str' 'nUnits_minus_str' 
    
     We also want to keep track of mappable length, in order to
     later compute mappability-adjusted RPKM for the whole feature type.
     Therefore gene_entry() also gets 'mapLength'.

     Finally, we want to keep track of gene boundaries. For this, we need
     keep track of the outermost ends of all features of a gene. Therefore
     gene_entry gets genomic coordinates of the region 'chrom' 'start' 'end'

    PROCEDURE OF READ DISTRIBUTION:
    Since we know 'geneStrand' strandedness of the gene, and total number of
    genes that overlap the region ('nUnits_plus_str' 'nUnits_minus_str'), we
    can correctly distribute reads between genes in the region.

    If a region contains overlapping features of genes transcribed from the
    same strand then reads aligned to plus and minus strands will be
    distributed evenly as counts of sense and antisense reads of the
    overlapping genes.

    If a region contains overlapping features of genes transcribe from the
    opposite genomic strands, the the reads aligned to plus and minus strands
    will be distributed to be sense reads of the two genes.

    '''

    if geneID in summary:
        pass
    else:
        summary = initialize_geneID( summary, geneID )
   
    # keep track of gene bounaries
    summary[geneID]['chrom']  = chrom
    summary[geneID]['strand'] = geneStrand
    if partType == 'exon' or partType == 'utr5p' or partType == 'utr3p':
        if summary[geneID]['start'] == '':
            summary[geneID]['start'] = start
        elif summary[geneID]['start'] > start:
            summary[geneID]['start'] = start
        else:
            pass
        if summary[geneID]['end'] == '':
            summary[geneID]['end'] = end
        elif summary[geneID]['end'] < end:
            summary[geneID]['end'] = end
        else:
            pass

    if geneStrand == '+' :
        summary[geneID]['sense'][partType]["MAPLENGTH"]     += mapLength 
        summary[geneID]['antisense'][partType]["MAPLENGTH"] += mapLength 
        summary[geneID]['sense'][partType]["COUNT"]         += \
                                                       pCount / nUnits_plus_str
        if nUnits_minus_str == 0:
            summary[geneID]['antisense'][partType]["COUNT"] += \
                                                       mCount / nUnits_plus_str
    elif geneStrand == '-':
        summary[geneID]['sense'][partType]["MAPLENGTH"]     += mapLength 
        summary[geneID]['antisense'][partType]["MAPLENGTH"] += mapLength 
        summary[geneID]['sense'][partType]["COUNT"]         += \
                                                      mCount / nUnits_minus_str
        if nUnits_plus_str == 0:
            summary[geneID]['antisense'][partType]["COUNT"] += \
                                                      pCount / nUnits_minus_str
    else:
        sys.exit("<%s> FATAL ERROR: could not parse strand info \
        !!!" % sys.argv)

    return summary


def count_of_units(units, repetetive_genes):
    '''
    Count number of genes (units) that are on forward and reverse strands.

    Note, if repetitive_genes are present, then only actual repetitive genes
    (i.e. those with '#' in their name) are counted.

    '''
    nUnits_plus_str  = 0
    nUnits_minus_str = 0
    for unit in units:
        myMatch  = re.match('(.*)~(.*){(.)}', unit)
        partType = myMatch.group(1)
        geneID   = myMatch.group(2)
        strand   = myMatch.group(3)
        if repetetive_genes == "present" and re.search('#', geneID):
            if strand == '+':
                nUnits_plus_str  += 1
            elif strand  == '-':
                nUnits_minus_str += 1
            else:
                sys.exit('<%s> FATAL ERROR: unable to interpret strand !!!' \
                % sys.argv[0])
        elif repetetive_genes == "absent":
            if strand == '+':
                nUnits_plus_str  += 1
            elif strand  == '-':
                nUnits_minus_str += 1
            else:
                sys.exit('<%s> FATAL ERROR: unable to interpret strand !!!' \
                % sys.argv[0])
    return (nUnits_plus_str, nUnits_minus_str)


def summarize_per_gene( plusTableFileObj, minusTableFileObj, mappedReads ):
    '''
    Read in two input files at the same time (this may be changed to reading
    one file at a time, as it will allow to work with regular RNA-seq, that is
    not strand specific). 
    
    On each, line count number of genes (units) overlapping the region with
    count_of_units() function. Counts of genes on (+) and (-) genomic strands
    stored separately in 'nUnits_plus_str' 'nUnits_minus_str'

    Check whether any non-coding and or repetitive genes overlap with the
    region by trying to match '# among the gene names. (the repetitive genes
    are distinguished by '#' in their name see
    "EXPRESSION_BY_FEATURE/bin/fix_genic_features.py"). If '#', set
    repetitive_genes = "present".
   
    Loop through each gene (unit) in the region and distribute counts between
    with gene_entry() function. 
    
    Note that if '#' is present (i.e. repetitive non-coding genes), then
    distributing of reads happens only among genes (units) with '#' in their
    name while regular genes are skipped.

    ''' 

    summary = {} 
    for plus_line in plusTableFileObj:
        
        minus_line =  minusTableFileObj.readline()
       
        pChrom, pStart, pEnd, pName, pRpkm, pCount, pLength, pMapLength \
                                                      = parse_line( plus_line )
        mChrom, mStart, mEnd, mName, mRpkm, mCount, mLength, mMapLength \
                                                     = parse_line( minus_line )
                                                                                                                                  
        assert pChrom  == mChrom, \
                       "<%s> FATAL ERROR: lines in -p <plusTab> -m <minusTab> \
                              are not matching by chromosome !!!" % sys.argv[0]
        chrom = pChrom                                                                                                            
        assert pStart  == mStart, \
                       "<%s> FATAL ERROR: lines in -p <plusTab> -m <minusTab> \
                                   are not matching by start !!!" % sys.argv[0]
        start = pStart                                                                                                            
        assert pEnd    == mEnd, \
                       "<%s> FATAL ERROR: lines in -p <plusTab> -m <minusTab> \
                                     are not matching by end !!!" % sys.argv[0]
        end   = pEnd                                                                                                              
        assert pName   == mName, \
                       "<%s> FATAL ERROR: lines in -p <plusTab> -m <minusTab> \
                                    are not matching by name !!!" % sys.argv[0]
        assert pMapLength == mMapLength, \
                       "<%s> FATAL ERROR: lines in -p <plusTab> -m <minusTab> \
                         are not matching by mappable length !!!" % sys.argv[0]
        mapLength = pMapLength

        # skip processing the row if annoation is 'intergenic'
        if pName == 'intergenic':
            continue
        
        # Check whether repetitive genes are present in annotation.
        # The repetitive genes are distinguished by '#' in their name.
        # See
        # "GENOME_TOOLS/annotate_genome_with_utrs_exons_introns.PerGene_SS.py"
        repetitive_genes = "absent" 
        if re.search('#', pName):
            repetitive_genes = "present"
        
        #
        # loop through all individual units/genes of name-annotation 
        #
        
        units = pName.split(";")

        (nUnits_plus_str, nUnits_minus_str) = \
                                        count_of_units(units, repetitive_genes)

        for unit in units:
            myMatch     = re.match('(.*)~(.*){(.)}', unit)
            partType    = myMatch.group(1)
            geneID      = myMatch.group(2)
            geneStrand  = myMatch.group(3)

            if repetitive_genes == "present":
                if re.search('#', geneID):
                    summary = gene_entry( summary, geneID, chrom, start, end, \
                                          geneStrand, partType, mapLength, \
                                          pCount, mCount, nUnits_plus_str, \
                                          nUnits_minus_str )
                else:
                    # when repetitive gene is present, don't give other genes
                    # any counts, but still pass to 'gene_entry()' in order
                    # to keep track of feature length and gene boundaries 
                    summary = gene_entry( summary, geneID, chrom, start, end, \
                                          geneStrand, partType, mapLength, \
                                          0, 0, 1, 1 )
            elif repetitive_genes == "absent":
                summary = gene_entry( summary, geneID, chrom, start, end, \
                                      geneStrand, partType, mapLength, \
                                      pCount, mCount, nUnits_plus_str, \
                                      nUnits_minus_str )
    return summary


def summarize_per_partition(plusTableFileObj, minusTableFileObj, mappedReads):
    summary = {}
    return summary


def summarize_table( plus_tabfile, minus_tabfile, tableType, mappedReads ):
    plusTableFileObj  = open(plus_tabfile,  mode = 'r')
    minusTableFileObj = open(minus_tabfile, mode = 'r')
    if tableType == "gene":
        summary = summarize_per_gene( plusTableFileObj, \
                                               minusTableFileObj, mappedReads )
    elif tableType == "partition":
        summary = summarize_per_partition( plusTableFileObj, \
                                               minusTableFileObj, mappedReads )
    else:
        sys.exit("FATAL ERROR: could not parse table type!!!")
    return summary


def print_header( partitions_list ):
    '''
    Print the header of the table
    '''
    outline = []
    outline.append('chrom\tstart\tend\tname\tannotation\tstrand')
    for partition in partitions_list:
        outline.append(partition + '_sense_mapRPKM')
        outline.append(partition + '_sense_reads')
        outline.append(partition + '_antisense_mapRPKM')
        outline.append(partition + '_antisense_reads')
        outline.append(partition + '_mapLen')
    outline.append('all_exons_sense_mapRPKM')
    outline.append('all_exons_sense_reads')
    outline.append('all_exons_sense_mapLen')
    print('\t'.join(outline), end = '\n') 


def get_rpkm_and_count(count, length, mappedReads):
    if length == 0:
        rpkm  = 'NA'
        count = 'NA'
    else:
        rpkm = (10**9 * count) / (length * mappedReads)
    return (rpkm, count)


def print_summary_rpkm( summary, mappedReads ):
    '''
    Add up counts of sense and antisense reads for features of the same type in
    each gene. 
    
    Add up mappabile lengths of features of the same type in each gene to
    compute mappability adjusted RPKM. 
    
    Add up counts of sense reads of 'utr5p', 'exon' and 'utr3p' features. Also,
    add up mappable lengths of these features. In the end produce one read
    count and one mappability adjusted RPKM statistic to characterize sense
    exonic expression of every gene.
    '''

    for geneID in sorted(summary):
        
        partitionVals = []
        partitionVals.append(summary[geneID]['chrom'])
        partitionVals.append(str(summary[geneID]['start']))
        partitionVals.append(str(summary[geneID]['end']))
        partitionVals.append(geneID)
        partitionVals.append(summary[geneID]['annotation'])
        partitionVals.append(summary[geneID]['strand'])
        all_sense_exons_reads  = 0
        all_sense_exons_mapLen = 0

        for partition in PARTITIONS_FOR_GENES:
            senseCount   = summary[geneID]['sense'][partition]["COUNT"]
            senseLength  = summary[geneID]['sense'][partition]["MAPLENGTH"]
            antiCount    = summary[geneID]['antisense'][partition]["COUNT"]
            antiLength   = summary[geneID]['antisense'][partition]["MAPLENGTH"]
            assert senseLength == antiLength, \
            "<%s> FATAL ERROR: lines in -p <plusTab> -m <minusTab> \
            are not matching by mappable length !!!" % sys.argv[0]
            if partition == 'exon' or partition == 'utr5p' \
                                                       or partition == 'utr3p':
                all_sense_exons_reads  += senseCount
                all_sense_exons_mapLen += senseLength
            
            (senseRpkm, senseCount) = \
                       get_rpkm_and_count(senseCount, senseLength, mappedReads)
            (antiRpkm,  antiCount)  = \
                         get_rpkm_and_count(antiCount, antiLength, mappedReads)

            partitionVals.append(str(senseRpkm))
            partitionVals.append(str(senseCount))
            partitionVals.append(str(antiRpkm))
            partitionVals.append(str(antiCount))
            partitionVals.append(str(senseLength))

        (all_sense_exons_mapRPKM, all_sense_exons_count) = \
            get_rpkm_and_count(all_sense_exons_reads, all_sense_exons_mapLen, \
                                                                   mappedReads)
        partitionVals.append(str(all_sense_exons_mapRPKM))
        partitionVals.append(str(all_sense_exons_count))
        partitionVals.append(str(all_sense_exons_mapLen))
        print('\t'.join(partitionVals), end = '\n')


def print_table( summary, mappedReads ):
    print_header( PARTITIONS_FOR_GENES )
    print_summary_rpkm( summary, mappedReads)



def main():
   
    plus_tab_file, minus_tab_file, mapped_reads, table_type = \
                                                         getOptions( sys.argv )
    mappedreads = parse_mappedReads( mapped_reads )
    
    summary = summarize_table( plus_tab_file, minus_tab_file, table_type, \
                                                                  mappedreads )
    print_table( summary, mappedreads )



if __name__ == '__main__':
    main()
