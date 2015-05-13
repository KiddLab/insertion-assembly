import genutils
import sys
import os.path
import numpy as np


from optparse import OptionParser

###############################################################################
def get_calls_from_vcf(vcfFileName):
    calls = []
    inFile = open(vcfFileName,'r')
    for line in inFile:
        if line[0] == '#':
           continue
        line = line.rstrip()
        line = line.split()
        c = line[0]
        p = int(line[1])
        #changed above for HGDP to accomodate the adding of begin/end coord's in get-filters-levels.py
        calls.append([c,p])
    inFile.close()
    return calls
###############################################################################
USAGE = """
python gather-soft-clipped-persample.py --vcf <input vcf file> --bam <input bam file> --out <output read file>

Goes through VCF and identifies softclip reads near breakpoint
Imposes minimum filtering critiera


"""
parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcfFile', help = 'input vcf file name')
parser.add_option('--bam',dest='bamFile', help = 'input bam file name')
parser.add_option('--out',dest='outFile', help = 'output file for soft clipped reads')


(options, args) = parser.parse_args()

if options.vcfFile is None:
    parser.error('vcf file name not given')
if options.bamFile is None:
    parser.error('BAM file name not given')
if options.outFile is None:
    parser.error('reads output file name not given')
###############################################################################


inputVCFCalls = options.vcfFile
inputBAMFile = options.bamFile
outPutReadFile = options.outFile

# hardcoded values
min_mean_qual = 20
minSoftClip = 20
window_size= 200
max_check_pair_size_for_soft_clip = 3000
min_map_q = 20
offSet = 33

# will just add readlen to the starts to get intervals, this will be off by 1, but will match
# the output of RetroSeq

# get the calls
calls = get_calls_from_vcf(inputVCFCalls)
print 'Read in %i calls from VCF' % (len(calls))

outFile = open(outPutReadFile,'w')
rN = 0
for call in calls:
    rN += 1
    if rN % 100 == 0:
        print 'Doing %i of %i ...' % (rN,len(calls))
    c = call[0]
    p = call[1]
    b = p - window_size
    e = p + window_size
    reg = '%s:%i-%i' % (c,b,e)
#    print reg
    numReads = 0
    numSC = 0
    bamIn = genutils.open_bam_read(inputBAMFile,reg)
    for line in bamIn:
        numReads += 1
        line = line.rstrip()
        line = line.split()
        samParse = genutils.parse_sam_line(line)

        #get rid of things that are not good
        if samParse['unMapped'] is True:
            continue
        if samParse['isDuplicate'] is True:
            continue
        if samParse['mapQ'] < min_map_q:
            continue
        if samParse['isPaired'] is False:
            continue
        if samParse['mateUnmapped'] is True:
            continue
        
        # check min soft clip
        if samParse['cigarCounts']['S'] < minSoftClip:
            continue
                
        qual = samParse['qual']
        cigar = samParse['cigar']
        cExpand = genutils.expand_cigar(cigar)
        if cExpand[0][1] == 'S':
            n = cExpand[0][0]
            qualPart = qual[0:n]        
        elif cExpand[-1][1] == 'S':
            n = cExpand[-1][0]
            qualPart = qual[(len(qual) - n) : ]
        else:
            print 'qhat??'
            print line
            print qual,cigar
            print cExpand
            sys.exit()

        #check that soft clipped portion passes min qual
        
        qualList = []
        for c in qualPart:
            i = ord(c) - offSet
            qualList.append(i)
        m = np.mean(qualList)
        if m < min_mean_qual:
            continue

        
        # check if ins size is reasonable and has min softClip values
        # not checking if softclip is near breakpoint, as will take all and hope that the
        # assembler will figure it out for us
        if (samParse['mateChrom'] == '=') and (abs(samParse['fragLen']) <= max_check_pair_size_for_soft_clip) and (samParse['cigarCounts']['S'] >= minSoftClip):
            #print line
            numSC += 1
            c = samParse['chrom']
            b = samParse['chromPos']
            e = b + samParse['seqLen']
            n = samParse['seqName']
            outFile.write('%s\t%i\t%i\tSOFTCLIP\t%s\t%s\n' % (c,b,e,n,samParse['cigar']))
    bamIn.close()
outFile.close()
