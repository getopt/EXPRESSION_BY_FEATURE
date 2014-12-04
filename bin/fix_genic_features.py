#!/usr/bin/env python

''' Usage:
    -t or --gtf        <file>           < transcript annotation in GTF format
    
    -g or --genomeFile <file>           < standard genome file with chromosome
                                          sizes
    
    -a or --assembly   <dm3|mm10|hg19>  < specify genome assembly in order to
                                          correctly gene names of non-coding
                                          RNAs

DESCRIPTION: 
Genome annotation for this list of genic features: 
    upstream|downstream|utr5p|utr3p|exon|intron

TODO:
    issue with per-gene dictionary and exon sorting
    clean up 'get_intron' by sorting the dictionary appropriately

'''

from __future__ import print_function
import pysam
import getopt
import sys
import re
import collections
import operator

PRIORITIES       =  \
         ['utr5p','utr3p','exon','intron','downstream','upstream','intergenic'] 
UPSTREAM_LEN     =  1000 # number of bp to consider to be the upstream region
DOWNSTREAM_LEN   =  1000 # number of bp to consider to be the upstream region
ASSEMBLY         =  ''

def getOptions(argv):
    '''
    Parse command line arguments
    '''
    gtfFile          = ''
    genomeFile       = ''
    assembly         = ''
    try:
        opts, args = getopt.getopt(argv[1:], "ht:g:a:", \
                                  ["help", "gtf=", "genomeFile=", "assembly="])
    except getopt.GetoptError:
        print('Usage: %s -t <gtf> -g <genomeFile> \
                    -a <assembly:dm3|mm10|hg19>' % argv[0], file = sys.stderr ) 
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('Usage: %s -g <gtf> -g <genomeFile> \
                    -a <assembly:dm3|mm10|hg19>' % argv[0], file = sys.stderr ) 
            sys.exit()
        elif opt in ("-t", "--gtf"):
            gtfFile  = arg
        elif opt in ("-g", "--genomeFile"):
            genomeFile  = arg
        elif opt in ("-a", "--assembly"):
            assembly  = arg
        else:
            assert False, "unhandled option"
    return gtfFile, genomeFile, assembly

def parse_genome_file ( genomefile ):
    '''
    Store chromosome names and lengths in a dictionary
    '''
    genomeFileObj = open(genomefile, mode = 'r')
    genomeDict = {}
    for line in genomeFileObj:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        size = int(fields[1])
        genomeDict[chrom] = size
    print( 'Finished reading geneome file', file =  sys.stderr ) 
    return genomeDict

def alter_if_noncoding(geneName, transcriptID, start, end):
    '''
    Identify names of "suspicious" (usually non-coding) genes (depends on
    assembly), and, when non-coding, joing 'gene_name' and 'transcript_id' on
    '#' attributes to:
    - later distinguish non-coding genes
    - make geneName trully unique (otherwise a problem for miRNA genes).
    
    NOTE: in dm3, matching on ':' also matches histone genes 
    '''
    if (ASSEMBLY == 'dm3'  and ( len(re.findall('mir-', geneName)) > 0  \
                              or len(re.findall('RNA:', geneName)) > 0 ) \
                              or len(re.findall(':', geneName)) > 0 \
                              or len(re.findall('^CR', geneName)) > 0 ) \
       or (ASSEMBLY == 'mm10' and ( len(re.findall('Mir',  geneName)) > 0 
                                 or len(re.findall('Snor', geneName))  > 0 )) \
       or (ASSEMBLY == 'hg19' and ( len(re.findall('SNOR', geneName)) > 0 \
                                 or len(re.findall('MIR',  geneName))  > 0 )):
        geneName = geneName + '#' + transcriptID
        # extend non-coding by 25nt on both sides for conservative counting
        start = start - 25 
        end   = end   + 25
    return (geneName, start, end)

def parse_gtf ( gtffile ):
    ''' [construct gtfDict]
    For every distict `gene_name` store all corresponding features from the GTF
    input file (i.e. all exons, and also all start and stop codons in case of
    coding genes).
    '''
    gtfObj = open(gtffile, mode = 'r')
    gtfDict = {}
    for line in gtfObj:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom         =  fields[0]
        feature       =  fields[2]
        start         =  int(fields[3]) - 1  # NOTE: GTF is always 1-based
        end           =  int(fields[4])
        strand        =  fields[6]
        
        # Special characters are temporarily removed from gene names to
        # ease use of gene names for pattern matching in 'get_strand' and 
        # 'get_featTypes_of_gene'
        fields[8]     =  re.sub('\(', '__paranthesisLeft__',    fields[8]) 
        fields[8]     =  re.sub('\)', '__paranthesisRight__',   fields[8]) 
        fields[8]     =  re.sub('\[', '__sqParanthesisLeft__',  fields[8]) 
        fields[8]     =  re.sub('\]', '__sqParanthesisRight__', fields[8]) 
        fields[8]     =  re.sub("\'", "__singleQuote__", fields[8])

        transcriptID  =  fields[8].split('transcript_id "')[1].split('"')[0]
        assert len(fields[8].split('gene_name "')) > 1, \
        "FATAL ERROR: gene_name field is absent from GTF input file !!!"
        geneName      =  fields[8].split('gene_name "')[1].split('"')[0]
        
        (geneName, start, end) = alter_if_noncoding(geneName, transcriptID,\
                                                                    start, end)

        if gtfDict.has_key(chrom):
            pass
        else:
            gtfDict[chrom] = {}
        if gtfDict[chrom].has_key(geneName):
            pass
        else:
            gtfDict[chrom][geneName] = {}
        if gtfDict[chrom][geneName].has_key(feature):
            gtfDict[chrom][geneName][feature].append(( chrom, start, end, \
                                                                      strand ))
        else:
            gtfDict[chrom][geneName][feature] = []
            gtfDict[chrom][geneName][feature].append(( chrom, start, end, \
                                                                      strand ))
    print( 'Finished first parsing GTF file', file = sys.stderr )  
    return gtfDict

def get_utr5p( gtfDict, chrom, gene, strand ):
    ''' [contribute to featureDict]
    Select a single *outermost* start-codon per gene as a reference position.
    Designate exons (or parts of exons) into 'utr5p' feature if they
    lay outside of the reference position.
    '''
    utr5p_list = []
    if gtfDict[chrom][gene].has_key('start_codon'):
        if strand == '+':
            start_codon_list = sorted(gtfDict[chrom][gene]['start_codon'])
            utr5p_end = start_codon_list[0][1] # i.e. left-most start codon
            myStart = ''
            myEnd   = ''
            for i in range(len(gtfDict[chrom][gene]['exon'])):
                if utr5p_end > gtfDict[chrom][gene]['exon'][i][1]:
                    myStart = gtfDict[chrom][gene]['exon'][i][1]
                    if utr5p_end >= gtfDict[chrom][gene]['exon'][i][2]:
                        myEnd = gtfDict[chrom][gene]['exon'][i][2]
                    else:
                        myEnd = utr5p_end
                    utr5p_list.append((myStart,myEnd,strand))
        elif strand == '-':
            start_codon_list = sorted(gtfDict[chrom][gene]['start_codon'])
            utr5p_begin = start_codon_list[-1][2] # i.e. right-most start codon
            myStart = ''
            myEnd   = ''
            for i in range(len(gtfDict[chrom][gene]['exon'])):
                if utr5p_begin < gtfDict[chrom][gene]['exon'][i][2]:
                    myEnd   = gtfDict[chrom][gene]['exon'][i][2]
                    if utr5p_begin >= gtfDict[chrom][gene]['exon'][i][1]:
                        myStart = utr5p_begin
                    else:
                        myStart = gtfDict[chrom][gene]['exon'][i][1]
                    utr5p_list.append((myStart,myEnd,strand))
    return utr5p_list

def get_utr3p( gtfDict, chrom, gene, strand ):
    ''' [contribute to featureDict]
    Select a single *outermost* stop-codon per gene as a reference position.
    Designate exons (or parts of exons) into 'utr3p' feature if they
    lay outside of the reference position.
    '''
    utr3p_list = []
    if gtfDict[chrom][gene].has_key('stop_codon'):
        if strand == '+':
            stop_codon_list = sorted(gtfDict[chrom][gene]['stop_codon'])
            utr3p_begin = stop_codon_list[-1][2] # i.e. right-most stopcodon
            myStart = ''
            myEnd   = ''
            for i in range(len(gtfDict[chrom][gene]['exon'])):
                if utr3p_begin < gtfDict[chrom][gene]['exon'][i][2]:
                    myEnd   = gtfDict[chrom][gene]['exon'][i][2]
                    if utr3p_begin >= gtfDict[chrom][gene]['exon'][i][1]:
                        myStart = utr3p_begin
                    else:
                        myStart = gtfDict[chrom][gene]['exon'][i][1]
                    utr3p_list.append((myStart,myEnd,strand))
        elif strand == '-':
            stop_codon_list = sorted(gtfDict[chrom][gene]['stop_codon'])
            utr3p_end = stop_codon_list[0][1] # i.e. left-most stopcodon
            myStart = ''
            myEnd   = ''
            for i in range(len(gtfDict[chrom][gene]['exon'])):
                if utr3p_end > gtfDict[chrom][gene]['exon'][i][1]:
                    myStart = gtfDict[chrom][gene]['exon'][i][1]
                    if utr3p_end >= gtfDict[chrom][gene]['exon'][i][2]:
                        myEnd = gtfDict[chrom][gene]['exon'][i][2]
                    else:
                        myEnd = utr3p_end
                    utr3p_list.append((myStart,myEnd,strand))
    return utr3p_list


def get_intron(gtfDict, chrom, gene, strand ):
    ''' [contribute to featureDict]
    Loop through every exon of a gene. Select regions between exons and
    designate them as 'intron'.
    NOTE: ISSUE WITH PER-GENE DICTIONARY AND EXON SORTING
    '''
    intron_list = []
    myStart = ''
    myEnd   = ''
    if len(gtfDict[chrom][gene]['exon']) == 1:
        return intron_list
    for i in range( 1, len(gtfDict[chrom][gene]['exon']) ):
        myStart = gtfDict[chrom][gene]['exon'][i - 1][2]
        myEnd   = gtfDict[chrom][gene]['exon'][i][1]
        intron_list.append((myStart,myEnd,strand))
    return intron_list

def get_exon(gtfDict, chrom, gene, strand ):
    ''' [contribute to featureDict]
    Easily transfer 'exon' features from gtfDict to featureDict.
    '''
    exon_list = []
    myStart = ''
    myEnd   = ''
    for i in range(len(gtfDict[chrom][gene]['exon'])):
        myStart = gtfDict[chrom][gene]['exon'][i][1]
        myEnd   = gtfDict[chrom][gene]['exon'][i][2]
        exon_list.append((myStart,myEnd,strand))
    return exon_list

def get_upstream(gtfDict, chrom, gene, strand, genomeDict):
    ''' [contribute to featureDict]
    Since you know all the exons, find the *outermost* and get 1kb upstream as
    the upstream region
    '''
    upstream_list = []
    myStart = ''
    myEnd   = ''
    if strand == '+':
        myEnd = sorted(gtfDict[chrom][gene]['exon'])[0][1]
        if myEnd - UPSTREAM_LEN > 0:
            myStart = myEnd - UPSTREAM_LEN
        else:
            myStart = 0
    elif strand == '-':
        nExons = len(gtfDict[chrom][gene]['exon'])
        myStart = sorted(gtfDict[chrom][gene]['exon'], \
                                   key = operator.itemgetter(2))[nExons - 1][2]
        if myStart + UPSTREAM_LEN > genomeDict[chrom]:
            myEnd = genomeDict[chrom]
        else:
            myEnd = myStart + UPSTREAM_LEN
    upstream_list.append((myStart,myEnd,strand))
    return upstream_list

def get_downstream(gtfDict, chrom, gene, strand, genomeDict):
    ''' [contribute to featureDict]
    Since you know all the exons, find the *outermost* and get 1kb downstream
    as the downstream region
    '''
    downstream_list = []
    myStart = ''
    myEnd   = ''
    if strand == '+':
        nExons = len(gtfDict[chrom][gene]['exon'])
        myStart = sorted(gtfDict[chrom][gene]['exon'], \
                                   key = operator.itemgetter(2))[nExons - 1][2]
        if myStart + UPSTREAM_LEN > genomeDict[chrom]:
            myEnd = genomeDict[chrom]
        else:
            myEnd = myStart + UPSTREAM_LEN
    elif strand == '-':
        myEnd = sorted(gtfDict[chrom][gene]['exon'])[0][1]
        if myEnd - UPSTREAM_LEN > 0:
            myStart = myEnd - UPSTREAM_LEN
        else:
            myStart = 0
    downstream_list.append((myStart,myEnd,strand))
    return downstream_list

def get_features( gtfDict, genomeDict ):
    ''' 
    Construct featureDict by calling functions to get UTRs, introns, etc. and
    pile them up on one dictionary ('chrom' and 'gene' as keys)
    '''
    featureDict = {}
    for chrom in sorted(gtfDict):
        print( 'Inputting gtf features for chromosome %s' \
                                                   % chrom, file = sys.stderr )
        featureDict[chrom] = {}
        for gene in gtfDict[chrom]:
            featureDict[chrom][gene] = {}
            strand = gtfDict[chrom][gene]['exon'][0][3]
            featureDict[chrom][gene]['utr5p']        =    \
                                      get_utr5p( gtfDict, chrom, gene, strand )
            featureDict[chrom][gene]['utr3p']        =    \
                                      get_utr3p( gtfDict, chrom, gene, strand )
            featureDict[chrom][gene]['intron']       =    \
                                     get_intron( gtfDict, chrom, gene, strand )
            featureDict[chrom][gene]['exon']         =    \
                                       get_exon( gtfDict, chrom, gene, strand )
            featureDict[chrom][gene]['upstream']     =    \
                       get_upstream( gtfDict, chrom, gene, strand, genomeDict )
            featureDict[chrom][gene]['downstream']   =    \
                     get_downstream( gtfDict, chrom, gene, strand, genomeDict )
    print( 'Finished parsing features', file = sys.stderr )  
    return featureDict

def add_labels_to_coverageDict( coverageDict, featureDict, chrom, priorities ):
    indexes = range(len(priorities))
    indexes.reverse()
    coverageDict[chrom] = {}
    if chrom not in featureDict:
        print( 'Found no exon information in GTF file for chromosome %s' % chrom, file = sys.stderr )
        return coverageDict   
    for gene in featureDict[chrom].keys():
        for f,val in enumerate(indexes):
            feature = priorities[val]
            if featureDict[chrom][gene].has_key(feature):  
                thisFeature = featureDict[chrom][gene][feature]
                for i in range(len(thisFeature)):
                    start  = thisFeature[i][0]
                    end    = thisFeature[i][1]
                    strand = thisFeature[i][2]
                    label  = feature + '~' + gene + '{' + strand + '}'
                    for j in range( start, end ):
                        if coverageDict[chrom].has_key(j):
                            coverageDict[chrom][j].append(label)
                        else:
                            coverageDict[chrom][j] = []
                            coverageDict[chrom][j].append(label)
    return coverageDict

def get_genes( labels ):
    genes = {} 
    for feature in labels:
        gene_strand = re.sub('..*~','',feature)  
        gene        = re.sub('{.*','',gene_strand)
        genes[gene] = 1
    return genes.keys()

def get_strand( gene, labels ):
    # for_search_only = gene
    # for_search_only = re.sub('\(', '\(', gene)
    # for_search_only = re.sub('\)', '\)', for_search_only)
    for feature in labels:
        #if re.search(for_search_only, feature):
        if re.search(gene, feature):
            strand = re.sub('.*{', '', feature)
            strand = re.sub('}',   '', strand)
            return strand
    
def get_featTypes_of_gene( gene, labels ):
    featTypes = {}
    # for_search_only = gene
    # for_search_only = re.sub('\(', '\(', gene)
    # for_search_only = re.sub('\)', '\)', for_search_only)
    for feature in labels:
        # if re.search(for_search_only, feature):
        if re.search(gene, feature):
            featType = re.sub('~.*','',feature)
            featTypes[featType] = 1
    return featTypes.keys()

def find_top_featureType( featTypes ):
    # PRIORITIES = ['utr5p','utr3p','exon','intron','downstream','upstream','intergenic']
    # print(featTypes)
    for featureDef in PRIORITIES:
        for featureType in featTypes:
            if re.search(featureDef, featureType):
                return featureType 

def parse_labels( labels, start, end, chrom ):
    # print(labels, end = " ::: ")
    # print((start, end, chrom), end = " ::: ")
    output   = []
    featDict = {}
    genes    = get_genes(labels)
    # print(genes)
    for gene in genes:
        if gene == 'intergenic':
            # print('intergenic')
            output.append('intergenic')
            pass
        else:
            # print(gene, labels)
            strand  = get_strand(gene, labels)
            featTypes  = get_featTypes_of_gene(gene, labels)
            topFeatType  = find_top_featureType(featTypes)
            # print(topFeatType)
            if featDict.has_key(topFeatType):
                pass
            else:
                featDict[topFeatType] = [] 
            # print(gene)
            featDict[topFeatType].append(gene + "{" + strand + "}")
    if featDict.has_key('utr5p') or featDict.has_key('utr3p') or featDict.has_key('exon') or featDict.has_key('intron'):
        if 'downstream' in featDict:
            del featDict['downstream']
        if 'upstream' in featDict:
            del featDict['upstream']
    for feature in sorted(featDict):
        for gene_and_strand in featDict[feature]:
            output.append(feature + '~' + gene_and_strand)
    currentLine = (chrom , str(start), str(end), ";".join(output))
    return currentLine

def fix_up_special_characters( label ):
    label =  re.sub('__paranthesisLeft__',    '(', label) 
    label =  re.sub('__paranthesisRight__',   ')', label) 
    label =  re.sub('__sqParanthesisLeft__',  '[', label) 
    label =  re.sub('__sqParanthesisRight__', ']', label) 
    label =  re.sub("__singleQuote__",       "\'", label)
    return label

def print_and_reset( outputLines, minLabel ):
    ## for debugging
    # print("New Call with minimal label %s" % minLabel)
    # for line in outputLines:
    #     print(line)
    # return []
    chrom    = outputLines[0][0]
    start    = outputLines[0][1]
    end      = outputLines[ len(outputLines) - 1 ][2]
    minLabel = outputLines[ len(outputLines) - 1 ][3]
    minLabel = fix_up_special_characters( minLabel )
    outline = chrom + '\t' + str(start) + '\t' + str(end) + '\t' + minLabel
    print(outline, end = '\n')
    return []

def print_coverageDict( coverageDict, genomeDict, priorities, chrom ):
    outputLines = []
    start  = 0

    if coverageDict[chrom].has_key(start):
        currentLine   = parse_labels(coverageDict[chrom][start], start, 0, chrom) 
        currentAnn    = coverageDict[chrom][start]
        minimal_label = currentLine[3] 
    else:
        currentAnn    = ['intergenic']
        minimal_label = ['intergenic']

    for i in range(genomeDict[chrom]):
        
        ### . . . 0 0 0 N
        if coverageDict[chrom].has_key(i) and currentAnn == ['intergenic']:
            currentLine = parse_labels(currentAnn, start, i, chrom) 

            outputLines.append( currentLine )
            outputLines = print_and_reset( outputLines, minimal_label )
            
            start = i
            currentAnn    = coverageDict[chrom][i]
            currentLine   = parse_labels(currentAnn, start, i, chrom)
            outputLines.append( currentLine )
            minimal_label = currentLine[3]

        ### . . . N N N M
        if coverageDict[chrom].has_key(i) and currentAnn != coverageDict[chrom][i]:
            currentLine = parse_labels(currentAnn, start, i, chrom) 
            outputLines.append( currentLine )
            
            currentLine = parse_labels(coverageDict[chrom][i], start, i, chrom) 
            if minimal_label == currentLine[3]:
                pass
            elif minimal_label != currentLine[3]:
                outputLines   = print_and_reset( outputLines, minimal_label  )
                start = i
                currentLine = parse_labels(coverageDict[chrom][i], start, i, chrom) 
                outputLines.append( currentLine )
                minimal_label = currentLine[3]
                currentAnn  = coverageDict[chrom][i]

        ### . . . N N N 0
        if currentAnn != ['intergenic'] and not coverageDict[chrom].has_key(i):
            currentLine = parse_labels(coverageDict[chrom][i - 1], start, i, chrom)
            outputLines.append( currentLine )
            outputLines = print_and_reset( outputLines, minimal_label )
            
            start = i
            currentAnn    = ['intergenic']
            currentLine   = parse_labels(currentAnn, start, i, chrom)
            outputLines.append( currentLine )
            minimal_label = currentLine[3]

        ### . . . | end of chrom
        if i == (genomeDict[chrom] - 1):
            currentLine = parse_labels(currentAnn, start, (i + 1), chrom)

            outputLines.append( currentLine )
            outputLines = print_and_reset( outputLines, minimal_label )

def print_coverage( genomeDict, featureDict, priorities ):
    coverageDict = {}
    for chrom in sorted(genomeDict):
        print( 'Working on chrom <%s>' % chrom, file = sys.stderr )
        print( 'Inputting exon information for chromosome %s' % chrom, file = sys.stderr )
        coverageDict = add_labels_to_coverageDict( coverageDict, featureDict, chrom, priorities )
        print( 'Printing annotation for chrom <%s>' % chrom, file = sys.stderr )
        print_coverageDict( coverageDict, genomeDict, priorities, chrom )
        print( 'Finished printing annotation for chromosome %s' % chrom, file = sys.stderr )
        try:
            del coverageDict[chrom]
        except KeyError:
            pass

def main(): 
    gtffile, genomefile, ASSEMBLY = getOptions( sys.argv )
    
    genomedict      =  parse_genome_file( genomefile )
    gtfdict         =  parse_gtf( gtffile )
    featuredict     =  get_features( gtfdict, genomedict )
    
    print_coverage( genomedict, featuredict, PRIORITIES )
    
    print( 'The program has finished successfully', file = sys.stderr )

if __name__ == '__main__':
    main()
