import genutils
import glob

from optparse import OptionParser
###############################################################################
USAGE = """
python combine-sam-and-evidence.py --vcf <input vcf file> --softclip <soft clip reads file>
                  --discovery <discovery reads file>  --sam <sam file of selected reads>
                  --outseq <output sequence file>  --outnames <output read names file>
                  --samplename <sample name for tracking>
                  

Gather reads from sample BAM file to be used in assembly.
Includes soft clip reads identified previously, as well as reads that directly support call
based on RetroSeq Discovery file


"""
parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcfFile', help = 'input vcf file name')
parser.add_option('--discover',dest='discFile',help = 'reads file from RetroSeq Discover stage')
parser.add_option('--softclip',dest='softClipFile',help = 'softclip read names file')
parser.add_option('--sam',dest='samFile', help = 'SAM file of all gathered reads')
parser.add_option('--outseq',dest='outSeqFile', help = 'file for output of read sequences')
parser.add_option('--outnames',dest='outNames', help = 'file for output of read names')
parser.add_option('--samplename',dest='sampleName', help = 'name of sample')




(options, args) = parser.parse_args()

if options.vcfFile is None:
    parser.error('vcf file name not given')
if options.discFile is None:
    parser.error('retroseq discovery reads file name not given')
if options.softClipFile is None:
    parser.error('soft clip read file name not given')
if options.samFile is None:
    parser.error('SAM file of all gathered reads not given')
if options.outSeqFile is None:
    parser.error('sequence output file not given')
if options.outNames is None:
    parser.error('read name output file not given')
if options.sampleName is None:
    parser.error('sample name not given')

###############################################################################

# read out sam, soft clip reads file, and discovery reads file
# prints out output sequence file name and output reads file name
# and does a windowed intersection

inputVCFCalls = options.vcfFile
samFile = options.samFile
softClipReadsFile = options.softClipFile
discoveryReadsFileCalls = options.discFile
outSeqName = options.outSeqFile
outReadsName = options.outNames
outReadsNameIntersect = outReadsName + '.w500'

sN = options.sampleName

window_size= 500 + 101 # to take any read intersect with softclip



outSeq = open(outSeqName,'w')
outReads = open(outReadsName,'w')



print 'Writing formated output read names to',outReadsName    
inFile = open(discoveryReadsFileCalls,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    line = line[0:5]
    line.append(sN)
    nl = '\t'.join(line) + '\n'
    outReads.write(nl)
inFile.close()

inFile = open(softClipReadsFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    line = line[0:5]
    line.append(sN)
    nl = '\t'.join(line) + '\n'
    outReads.write(nl)
inFile.close()

print 'Writing formated output read sequences to',outSeqName    

inFile = open(samFile,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    samRec = genutils.parse_sam_line(line)

    seq = samRec['seq']
    qual = samRec['qual']
    if samRec['reverseStrand'] is True:
        qual = qual[::-1]
        seq = genutils.revcomp(seq)
    if samRec['isFirst'] is True:
        i = 0
    else:
        i = 1
    nl = [samRec['seqName'],sN,str(i),seq,qual]
    nl = '\t'.join(nl) + '\n'
    outSeq.write(nl)
inFile.close()

outSeq.close()
outReads.close()

print 'Now doing window bed to get reads associated with each call'
print 'Output to',outReadsNameIntersect

cmd = 'windowBed -w %i -a %s -b %s > %s ' % (window_size,outReadsName,inputVCFCalls,outReadsNameIntersect)
print 'Getting discovery support calls within %i' % (window_size)
print cmd
genutils.runCMD(cmd)
print 'COMPLETE'




