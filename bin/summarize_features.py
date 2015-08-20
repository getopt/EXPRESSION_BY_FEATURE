#!/usr/bin/env python

''' Usage:

-b or --bothTab          <tab>  path/to/table_for_(+ & -)_genomic_strand
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


INPUT: one or more tables 
       
       in '*.per_geneSS_[plus|minus].tab' format (for reads mapped to the
       forward [plus] genomic strand, and one for reads mapped to the reverse
       [minus] genomic strand) or '*.per_geneNS_[plus|minus].tab' when library
       is unstranded. All libraries look like:
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
    bothTabfiles   = ''
    plusTabfiles   = ''
    minusTabfiles  = ''
    mappedReads   = ''
    tableType     = 'gene'
    try:
        opts, args = getopt.getopt(argv[1:], "hb:p:m:n:t:", \
        [ "help", "bothTab=", "plusTab=", "minusTab=", "mappedReads=", "type=" ])
    except getopt.GetoptError:
        print ('Usage: %s -b <bothTab> -p <plusTab> -m <minusTab> -n <mappedReads> \
                                      -t <type>'  % argv[0], file = sys.stderr) 
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print ('Usage: %s -b <bothTab> -p <plusTab> -m <minusTab> -n <mappedReads> \
                                       -t <type>' % argv[0], file = sys.stderr)
            sys.exit()
        elif opt in ("-b", "--bothTab"):
            bothTabfiles = arg
        elif opt in ("-p", "--plusTab"):
            plusTabfiles = arg
        elif opt in ("-m", "--minusTab"):
            minusTabfiles = arg
        elif opt in ("-n", "--mappedReads"):
            mappedReads = arg
        elif opt in ("-t", "--type"):
            tableType = arg
        else:
            assert False, "unhandled option"
    return bothTabfiles, plusTabfiles, minusTabfiles, mappedReads, tableType


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
    summary[geneID]['undef'] = {}
    summary[geneID]['undef']['upstream']                    = {}
    summary[geneID]['undef']['upstream']["MAPLENGTH"]       = 0
    summary[geneID]['undef']['upstream']["COUNT"]           = 0
    summary[geneID]['undef']['utr5p']                       = {}
    summary[geneID]['undef']['utr5p']["MAPLENGTH"]          = 0
    summary[geneID]['undef']['utr5p']["COUNT"]              = 0
    summary[geneID]['undef']['exon']                        = {}
    summary[geneID]['undef']['exon']["MAPLENGTH"]           = 0
    summary[geneID]['undef']['exon']["COUNT"]               = 0
    summary[geneID]['undef']['intron']                      = {}
    summary[geneID]['undef']['intron']["MAPLENGTH"]         = 0
    summary[geneID]['undef']['intron']["COUNT"]             = 0
    summary[geneID]['undef']['utr3p']                       = {}
    summary[geneID]['undef']['utr3p']["MAPLENGTH"]          = 0
    summary[geneID]['undef']['utr3p']["COUNT"]              = 0
    summary[geneID]['undef']['downstream']                  = {}
    summary[geneID]['undef']['downstream']["MAPLENGTH"]     = 0
    summary[geneID]['undef']['downstream']["COUNT"]         = 0

    summary[geneID]['annotation'] = ''
    summary[geneID]['annotation'] = interprete_geneID( geneID )
    
    summary[geneID]['action'] = {} 

    return summary


def interprete_geneID( geneID ):
    annotation = 'NA'
    (geneName, annotation) = geneID.split("#")
    return annotation


def gene_entry( summary, geneID, chrom, start, end, geneStrand, partType, mapLength, count, genomicStrand, nUnits_plus_str, nUnits_minus_str ):
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
    
    # based on 'genomicStrand' deside what type of count we are dealing with 
    pCount = 0
    mCount = 0
    bCount = 0
    if genomicStrand == 'plus': 
        pCount = count 
    elif genomicStrand == 'minus': 
        mCount = count 
    elif genomicStrand == 'both': 
        bCount = count 
    else:
        sys.exit("<%s> FATAL ERROR: could not parse library direction \
        !!!" % sys.argv)

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
    
    summary[geneID]['sense'][partType]["MAPLENGTH"]     += mapLength 
    summary[geneID]['antisense'][partType]["MAPLENGTH"] += mapLength 
    summary[geneID]['undef'][partType]["MAPLENGTH"]     += mapLength 

    if bCount > 0:
        summary[geneID]['undef'][partType]["COUNT"]         += \
                                  bCount / (nUnits_plus_str + nUnits_minus_str)

    if geneStrand == '+' :
        summary[geneID]['sense'][partType]["COUNT"]         += \
                                                       pCount / nUnits_plus_str
        if nUnits_minus_str == 0:
            summary[geneID]['antisense'][partType]["COUNT"] += \
                                                       mCount / nUnits_plus_str
    elif geneStrand == '-':
        summary[geneID]['sense'][partType]["COUNT"]         += \
                                                      mCount / nUnits_minus_str
        if nUnits_plus_str == 0:
            summary[geneID]['antisense'][partType]["COUNT"] += \
                                                      pCount / nUnits_minus_str
    else:
        sys.exit("<%s> FATAL ERROR: could not parse strand info!!!" % sys.argv)

    return summary


def count_of_units(units):
    '''
    Count number of genes (units) that are on forward and reverse strands.

    Note, if 'AR' repetitive_genes are present, then only actual repetitive
    genes (i.e. those with '#AR_' in their name) are counted.

    '''
    n_AR_Units_plus_str  = 0
    n_AR_Units_minus_str = 0
    
    n_other_nc_Units_plus_str  = 0
    n_other_nc_Units_minus_str = 0
    
    n_protein_coding_Units_plus_str  = 0
    n_protein_coding_Units_minus_str = 0

    for unit in units:
        myMatch  = re.match('(.*)~(.*){(.)}', unit)
        partType = myMatch.group(1)
        geneID   = myMatch.group(2)
        strand   = myMatch.group(3)
        if  re.search('#AR_', geneID):
            if strand == '+':
                n_AR_Units_plus_str  += 1
            elif strand  == '-':
                n_AR_Units_minus_str += 1
            else:
                sys.exit('<%s> FATAL ERROR: unable to interpret strand !!!' \
                % sys.argv[0])
        elif re.search('ncRNA_', geneID):
            if strand == '+':
                n_other_nc_Units_plus_str += 1
            elif strand  == '-':
                n_other_nc_Units_minus_str += 1
            else:
                sys.exit('<%s> FATAL ERROR: unable to interpret strand !!!' \
                % sys.argv[0])
        elif re.search('#protein_coding', geneID):
            if strand == '+':
                n_protein_coding_Units_plus_str += 1
            elif strand  == '-':
                n_protein_coding_Units_minus_str += 1
            else:
                sys.exit('<%s> FATAL ERROR: unable to interpret strand !!!' \
                % sys.argv[0])
        else:
            sys.exit('<%s> FATAL ERROR: do not understand type of unit!!!' \
            % sys.argv[0])
    return ( n_AR_Units_plus_str, n_AR_Units_minus_str, \
             n_other_nc_Units_plus_str, n_other_nc_Units_minus_str, \
             n_protein_coding_Units_plus_str, n_protein_coding_Units_minus_str )


def summarize_per_gene( tableFile, genomicStrand, mappedReads, summary ):
    '''
    Read in one input file at a time keeping track of genomic strand to which
    reads are aligned: 'plus' for forward strand, 'minus' for reverse strand,
    'both' when library is not strand specific and strand information is
    meaningless.
    
    The input is basically a list of intervals and read counts for the interval.

    By parsing the 'name' attribute of the interval in the input file, for
    each interval count number of genes (units) overlapping the region with
    count_of_units() function. Counts of genes on (+) and (-) genomic strands
    stored separately in 'nUnits_plus_str' 'nUnits_minus_str'

    Check whether any 'AR_*' or 'lncRNA_encoding_wkncRNA' or 'ncRNA_other'
    genes overlap with the region. If 'AR_*' is present (i.e. well-known
    non-coding genes), then distributing of reads happens only among the 'AR_*'
    genes (units) with 'AR_' in their name while regular genes are skipped
    (i.e. the 'AR_*' genes absorb all the reads from coding genes). 
    
    If overlap with 'lncRNA_encoding_wkncRNA' or 'ncRNA_other' is detected,
    then they won't get any reads as long as they overlap with a protein coding
    gene on the same strand.
    
    In case of overlapping genes are on transcribed from opposite strand, then
    reads will be distributed according to the strandedness.

    Loop through each gene (unit) in the region and distribute counts between
    with gene_entry() function. 
    
    ''' 
    
    tableFileObj = open(tableFile, mode = 'r')
    

    for line in tableFileObj:
        
        chrom, start, end, name, rpkm, count, length, mapLength \
                                                      = parse_line( line )

        # skip processing the row if annoation is 'intergenic'
        if name == 'intergenic':
            continue
        
        # Check what kind of repetitive genes is present in annotation.  The
        # repetitive genes that are distinguished by 'AR' in their name will
        # absorb the reads, "other" non-coding genes (i.e.
        # 'lncRNA_encoding_wkncRNA' or 'ncRNA_other') will not get any
        # reads if they overlap with a 'protein_coding' gene on the same strand.
        repetitive_genes = "" 
        if re.search('#AR_', name):
            repetitive_genes = "AR"
        elif re.search('ncRNA_', name):
            repetitive_genes = "other_nc"
        else:
            repetitive_genes = "absent"

        #
        # loop through all individual units/genes of name-annotation 
        #
        
        units = name.split(";")

        ( n_AR_Units_plus_str, n_AR_Units_minus_str, \
        n_other_nc_Units_plus_str, n_other_nc_Units_minus_str, \
        n_protein_coding_Units_plus_str, n_protein_coding_Units_minus_str ) = \
                                                                count_of_units(units)
        action = '.' 
        for unit in units:
            myMatch     = re.match('(.*)~(.*){(.)}', unit)
            partType    = myMatch.group(1)
            
            geneID      = myMatch.group(2)
            geneStrand  = myMatch.group(3)

            # this is to prevent genes on different chrms from having the same names 
            # "|chrom" is subsequently stripped during printing
            geneID      = geneID + "|" + chrom  + "|" + geneStrand 
                                                
            if repetitive_genes == "AR":
                if re.search('#AR_', geneID):
                    if count > 0 and ( (n_protein_coding_Units_plus_str + n_protein_coding_Units_minus_str) > 0 \
                                  or   (n_other_nc_Units_plus_str + n_other_nc_Units_minus_str) > 0 ):
                        action = 'absorption'
                    elif count > 0 and ( (n_AR_Units_plus_str + n_AR_Units_minus_str) > 1 ):
                        action = 'sharing'
                    else:
                        action = 'singleton'
                    summary = gene_entry( summary, geneID, chrom, start, end, \
                                          geneStrand, partType, mapLength, \
                                          count, genomicStrand, \
                                          n_AR_Units_plus_str, \
                                          n_AR_Units_minus_str )
                else:
                    if count > 0 and ( (n_protein_coding_Units_plus_str + n_protein_coding_Units_minus_str) > 0 \
                                  or   (n_other_nc_Units_plus_str + n_other_nc_Units_minus_str) > 0 ):
                        action = 'absorption'
                    summary = gene_entry( summary, geneID, chrom, start, end, \
                                          geneStrand, partType, mapLength, \
                                          0, genomicStrand, 1, 1 )

            elif repetitive_genes == "other_nc":
                if re.search('ncRNA_', geneID):
                    # scenario when ncRNA 
                    #   - there are sense reads
                    #   - there is no coding gene on the same strand
                    #   - ncRNA gets sense reads
                    if (geneStrand == '+' and genomicStrand == 'plus'  \
                                                         and n_protein_coding_Units_plus_str == 0) \
                    or (geneStrand == '-' and genomicStrand == 'minus' \
                                                            and n_protein_coding_Units_minus_str == 0):
                        if (n_other_nc_Units_plus_str + n_other_nc_Units_minus_str ) > 1: 
                            action = 'sharing'
                        else:
                            action = 'singleton'
                        summary = gene_entry( summary, geneID, chrom, start, end, \
                                              geneStrand, partType, mapLength, \
                                              count, genomicStrand, \
                                              n_other_nc_Units_plus_str, \
                                              n_other_nc_Units_minus_str)
                    else:
                        if (count > 0 and (n_protein_coding_Units_plus_str + n_protein_coding_Units_minus_str ) > 0):
                          action = 'capture'
                          summary = gene_entry( summary, geneID, chrom, start, end, \
                                              geneStrand, partType, mapLength, \
                                              0, genomicStrand, 1, 1 )
                        else:
                        # scenario when ncRNA 
                        #   - there are antisense reads
                        #   - there is no coding gene on the same (antisense) strand
                        #   - ncRNA gets antisense reads
                          if (n_other_nc_Units_plus_str + n_other_nc_Units_minus_str ) > 1:
                            action = 'sharing'
                          else:
                            action = 'singleton'
                          summary = gene_entry( summary, geneID, chrom, start, end, \
                                              geneStrand, partType, mapLength, \
                                              count, genomicStrand, \
                                              n_other_nc_Units_plus_str, \
                                              n_other_nc_Units_minus_str )

                elif re.search('#protein_coding', geneID):
                    if (geneStrand == '+' and genomicStrand == 'minus'  \
                                                             and n_other_nc_Units_minus_str > 0) \
                     or (geneStrand == '-' and genomicStrand == 'plus'  \
                                                             and n_other_nc_Units_plus_str > 0):
                        summary = gene_entry( summary, geneID, chrom, start, end, \
                                              geneStrand, partType, mapLength, \
                                              0, genomicStrand, 1, 1 )
                    else:
                        summary = gene_entry( summary, geneID, chrom, start, end, \
                                              geneStrand, partType, mapLength, \
                                              count, genomicStrand, \
                                              n_protein_coding_Units_plus_str, \
                                              n_protein_coding_Units_minus_str )
                        if count > 0: action = 'capture'
                else:
                    sys.exit('<%s> FATAL ERROR: unable to geneID with ncRNA !!!' \
                    % sys.argv[0])

            elif repetitive_genes == "absent":
                if (n_protein_coding_Units_plus_str + n_protein_coding_Units_minus_str) > 1:
                    if count > 0: action = 'sharing'
                else:
                    if count > 0: action = 'singleton'
                summary = gene_entry( summary, geneID, chrom, start, end, \
                                      geneStrand, partType, mapLength, \
                                      count, genomicStrand, \
                                      n_protein_coding_Units_plus_str, \
                                      n_protein_coding_Units_minus_str )

            else:
                sys.exit('<%s> FATAL ERROR: unable to geneID at all !!!' \
                % sys.argv[0])

            if action in summary[geneID]['action']:
                pass
            else:
                summary[geneID]['action'][action] = 1
                action = '.'
                summary[geneID]['action'][action] = 1

    return summary


def summarize_per_partition(tableFile, genomicStrand, mappedReads, summary):
    summary = {}
    return summary


def summarize_table(tableFile, genomicStrand, mappedReads, tableType, summary):
    if tableType == "gene":
        summary = summarize_per_gene( tableFile, \
                                         genomicStrand, mappedReads, summary )
    elif tableType == "partition":
        summary = summarize_per_partition( tableFile, \
                                         genomicStrand, mappedReads, summary )
    else:
        sys.exit("FATAL ERROR: could not parse table type!!!")
    return summary


def run_summarize( both_tabfiles, plus_tabfiles, minus_tabfiles, tableType, mappedReads ):
    
    summary = {}
    nFiles  = 0

    bothFiles  = both_tabfiles.split(',') 
    plusFiles  = plus_tabfiles.split(',') 
    minusFiles = minus_tabfiles.split(',') 
    

    for bothFile in bothFiles:
        if bothFile != '':
            nFiles += 1
            summary = summarize_table( bothFile, \
                                      'both', mappedReads, tableType, summary )

    for plusFile in plusFiles:
        if plusFile != '':
            nFiles += 1
            summary = summarize_table( plusFile, \
                                      'plus', mappedReads, tableType, summary )

    for minusFile in minusFiles:
        if minusFile != '':
            nFiles += 1
            summary = summarize_table( minusFile, \
                                      'minus', mappedReads, tableType, summary )
    
    return (summary, nFiles)


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
        outline.append(partition + '_undef_mapRPKM')
        outline.append(partition + '_undef_reads')
        outline.append(partition + '_mapLen')
    outline.append('all_exons_sense_mapRPKM')
    outline.append('all_exons_sense_reads')
    outline.append('all_exons_sense_mapLen')
    outline.append('actions')
    print('\t'.join(outline), end = '\n') 


def get_rpkm_and_count(count, length, mappedReads):
    if length == 0:
        rpkm  = 'NA'
        count = 'NA'
    else:
        rpkm = (10**9 * count) / (length * mappedReads)
    return (rpkm, count)


def print_summary_rpkm( summary, mappedReads, nFiles ):
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
        partitionVals.append(geneID.split('#')[0])

        # This is to strip chromosome from annotation.
        # Attaching chrm was necessary in 'gene_entry()'
        # in order to prevent genes on different chrms 
        # from having the same name.
        myMatch    = re.match('(.*)\|(.*)\|(.*)', summary[geneID]['annotation']) 
        annotation = myMatch.group(1)                                             
        partitionVals.append(annotation)                                          
                                                                                  
        partitionVals.append(summary[geneID]['strand'])
        all_sense_exons_reads  = 0
        all_sense_exons_mapLen = 0

        for partition in PARTITIONS_FOR_GENES:
            senseCount   = summary[geneID]['sense'][partition]["COUNT"]
            senseLength  = summary[geneID]['sense'][partition]["MAPLENGTH"] / nFiles
            antiCount    = summary[geneID]['antisense'][partition]["COUNT"]
            antiLength   = summary[geneID]['antisense'][partition]["MAPLENGTH"] / nFiles
            bothCount    = summary[geneID]['undef'][partition]["COUNT"]
            bothLength   = summary[geneID]['undef'][partition]["MAPLENGTH"] / nFiles
            
            if partition == 'exon' or partition == 'utr5p' \
                                                       or partition == 'utr3p':
                all_sense_exons_reads  += senseCount + bothCount
                all_sense_exons_mapLen += senseLength
            
            (senseRpkm, senseCount) = \
                       get_rpkm_and_count(senseCount, senseLength, mappedReads)
            (antiRpkm,  antiCount)  = \
                         get_rpkm_and_count(antiCount, antiLength, mappedReads)
            (bothRpkm,  bothCount)  = \
                         get_rpkm_and_count(bothCount, bothLength, mappedReads)

            partitionVals.append(str(senseRpkm))
            partitionVals.append(str(senseCount))
            partitionVals.append(str(antiRpkm))
            partitionVals.append(str(antiCount))
            partitionVals.append(str(bothRpkm))
            partitionVals.append(str(bothCount))
            partitionVals.append(str(senseLength))

        (all_sense_exons_mapRPKM, all_sense_exons_count) = \
            get_rpkm_and_count( all_sense_exons_reads, all_sense_exons_mapLen, \
                                                                    mappedReads)
        partitionVals.append(str(all_sense_exons_mapRPKM))
        partitionVals.append(str(all_sense_exons_count))
        partitionVals.append(str(all_sense_exons_mapLen))
        
        actions = [] 
        for action in sorted(summary[geneID]['action']):
            actions.append(action)
        partitionVals.append(str(','.join(actions)))

        print('\t'.join(partitionVals), end = '\n')


def print_table( summary, mappedReads, nFiles ):
    print_header( PARTITIONS_FOR_GENES )
    print_summary_rpkm( summary, mappedReads, nFiles)



def main():
   
    both_tab_files, plus_tab_files, minus_tab_files, mapped_reads, table_type = \
                                                         getOptions( sys.argv )
    mappedreads = parse_mappedReads( mapped_reads )
    
    (summary, nFiles) = run_summarize( both_tab_files, plus_tab_files, minus_tab_files, table_type, \
                                                                  mappedreads )
    
    print_table( summary, mappedreads, nFiles )



if __name__ == '__main__':
    main()
