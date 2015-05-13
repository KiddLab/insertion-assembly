import genutils
import sys
import os.path

from optparse import OptionParser
###############################################################################
USAGE = """
python gather-seq-sam-persample.py --vcf <input vcf file> --bam <input bam file>  --softclip <soft clip reads file>
                  --discovery <discovery reads file>  --out <output reads file>
                  

Gather reads from sample BAM file to be used in assembly.
Includes soft clip reads identified previously, as well as reads that directly support call
based on RetroSeq Discovery file


"""
parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcfFile', help = 'input vcf file name')
parser.add_option('--bam',dest='bamFile', help = 'input bam file name')
parser.add_option('--out',dest='outFile', help = 'output SAM file of gathered reads')
parser.add_option('--discover',dest='discFile',help = 'reads file from RetroSeq Discover stage')
parser.add_option('--softclip',dest='softClipFile',help = 'softclip read names file')


(options, args) = parser.parse_args()

if options.vcfFile is None:
    parser.error('vcf file name not given')
if options.bamFile is None:
    parser.error('BAM file name not given')
if options.discFile is None:
    parser.error('retroseq discovery reads file name not given')
if options.softClipFile is None:
    parser.error('soft clip read file name not given')
if options.outFile is None:
    parser.error('reads output SAM file name not given')
###############################################################################


# set up file info to get


inputBAMFile = options.bamFile
if os.path.isfile(inputBAMFile) is False:
    print 'FILE ERROR'
    print inputBAMFile
    sys.exit()


#inputVCFCalls = '/home/jwilds/kidd-lab/jwilds-projects/RetroSeq/HGDP/merged/HGDP.merged.HERV.redo.genome.notRef500.level8.vcf'
inputVCFCalls = options.vcfFile
softClipReadsFile = options.softClipFile
discoveryReadsFile = options.discFile
outputSAMFile = options.outFile
# file name of reads that support each call
discoveryReadsFileCalls = outputSAMFile + '.DICOVERY_SUPPORT'

# hard coded in, based on library insert size
window_size= 500



# get the discovery reads that actually support calls we are making

cmd = 'windowBed -w %i -a %s -b %s > %s ' % (window_size,discoveryReadsFile,inputVCFCalls,discoveryReadsFileCalls)
print 'Getting discovery support calls within %i' % (window_size)
print cmd
genutils.runCMD(cmd)

#read in reads to do
reads = {}
inFile = open(discoveryReadsFileCalls,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    rN = line[4]
    reads[rN] = 1
inFile.close()

print 'Have names of %i reads to gather, after discovery file ' % (len(reads))


inFile = open(softClipReadsFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    rN = line[4]
    reads[rN] = 1
inFile.close()

print 'Have names of %i reads to gather, after softclip file ' % (len(reads))



n = 0
outSAM = open(outputSAMFile,'w')
inBAM = genutils.open_bam_read(inputBAMFile)
for line in inBAM:
    ol = line
    line = line.split('\t')
    if line[0] in reads:
        outSAM.write(ol)
        n += 1
inBAM.close()
outSAM.close()

print 'Wrote out %i lines' % n
