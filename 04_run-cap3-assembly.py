import genutils
import os
import sys
from optparse import  OptionParser
###############################################################################
def make_fastq_file(readsInInterval,sN,nameToSeq,fastqFileName):
    outFile = open(fastqFileName,'w')
    for name in readsInInterval.keys():
        for i in range(2):
            outFile.write('@%s\n%s\n+\n%s\n' % (name,nameToSeq[(name,sN)][i][0],nameToSeq[(name,sN)][i][1]))    
    outFile.close()
###############################################################################
def run_repeatmasker(fastaFile):
    cmd = 'RepeatMasker --species human %s ' % (fastaFile)
    # change this to use other species or libraries
#    cmd = 'RepeatMasker --species dog %s ' % (fastaFile)
    if  os.path.isfile(fastaFile) is True:
        genutils.runCMD(cmd)
###############################################################################



#############################################################################

USAGE = """run-cap3-assembly.py --vcf <input vcf file> --seq <sequence file>  --namesintersect <names intersection file>
                                --outdirbase <base directory for assembly output>
                                --log  <log file name> 


"""
parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcfFile', help = 'input vcf file name')
parser.add_option('--seq',dest='seqFile', help = 'file formatted sequence reads')
parser.add_option('--namesintersect',dest='nameInter', help = 'file of names intersected with calls')
parser.add_option('--outdirbase',dest='outDirBase', help = 'base directory for assembly output')
parser.add_option('--log',dest='logFileName', help = 'log file name')


(options,args)=parser.parse_args()

if options.vcfFile is None:
    parser.error('vcf file name not given')
if options.seqFile is None:
    parser.error('sequence  file not given')
if options.nameInter is None:
    parser.error('read name intersect file not given')
if options.outDirBase is None:
    parser.error('base dir for assembly output not given')
if options.logFileName is None:
    parser.error('log file name not given')


    
#############################################################################


inputVCFCalls = options.vcfFile
seqFileName = options.seqFile
supportReadsFileIntersect = options.nameInter
logName = options.logFileName
outDirBase = options.outDirBase
window_size= 200

LOG = open(logName,'w',1)
print 'opened logname %s' % logName


if  os.path.isdir(outDirBase) is False:
    print 'Making dir...'
    cmd = 'mkdir -p ' + outDirBase
    print cmd
    genutils.runCMD(cmd)
else:
    print 'output dir exists'
    print outDirBase


if outDirBase[-1] != '/':
    outDirBase += '/'


intervalsToReadNames = {}
nameToSeq = {}
inFile = open(supportReadsFileIntersect,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    rN = line[4]
    sN = line[5]
    c = line[6]
    p = line[7]
    k = c + '_' + p
    if k not in intervalsToReadNames:
        intervalsToReadNames[k] = {}
    # seq,qual for rn1/rn2
    nameToSeq[(rN,sN)] = [['',''],['','']]
    if sN not in intervalsToReadNames[k]:
        intervalsToReadNames[k][sN] = {}
    intervalsToReadNames[k][sN][rN] = 1
    
inFile.close()
print 'Read in info for %i intervals and %i total fragments' % (len(intervalsToReadNames),len(nameToSeq))

inFile = open(seqFileName,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    rN = line[0]
    sN = line[1]
    rNum = int(line[2])
    seq = line[3]
    qual = line[4]
    if (rN,sN) in nameToSeq:
        nameToSeq[(rN,sN)][rNum] = [seq,qual]
inFile.close()
print 'read in sequences'



regionsToDo = []

inFile = open(inputVCFCalls,'r')
for line in inFile:
    if line[0] == '#':
        continue
    line = line.rstrip()
    line = line.split()
    c = line[0]
    p = line[1]
    k = c + '_' + p
    regionsToDo.append(k)
    
    #setup chrom
    outDirChrom = outDirBase + c
    if  os.path.isdir(outDirChrom) is False:
        print 'Making dir...'
        cmd = 'mkdir -p ' + outDirChrom
        print cmd
        genutils.runCMD(cmd)

inFile.close()


print 'Read in %i regions to do' % (len(regionsToDo))

for region in regionsToDo:
    print '***',region,'***'
    if region not in intervalsToReadNames:
        print 'FAIL'
        LOG.write('FAIL\t%s\n' % region)
        continue
    samplesInInterval = intervalsToReadNames[region]
    print 'Have %i samples to assemble' % len(samplesInInterval)
    parts = region.split('_')
    parts  = parts[0:-1]
    chrom = '_'.join(parts)

    outPutDir = outDirBase + chrom + '/' + region 
    if  os.path.isdir(outPutDir) is False:
        print 'Making dir...'
        cmd = 'mkdir -p ' + outPutDir
        print cmd
        genutils.runCMD(cmd)
    else:
        print 'output dir exists'
        print outPutDir

    outPutDir += '/'
    outReadsFile = outPutDir + 'reads.fq'    
    print outReadsFile
    fastqFile = open(outReadsFile,'w')

    for sN in samplesInInterval.keys():
        readsForSample =  samplesInInterval[sN] 
        numR = len(readsForSample.keys())
        print sN,numR
        for name in readsForSample.keys():
            for i in range(2):                
                #convert names
                capName = sN + '_' + name
                # Reads from SRA have a '.' in the name that will cause problems with CAP3
                capName = capName.replace('.','_')
                
                # for reads from SRA, add .1 and .2 read name
                if i == 0:
                    capName = capName + '.1'
                if i == 1:
                    capName = capName + '.2'               
                fastqFile.write('@%s\n%s\n+\n%s\n' % (capName,nameToSeq[(name,sN)][i][0],nameToSeq[(name,sN)][i][1]))    
    fastqFile.close()

    # convert to fasta
    cmd = 'fastq_to_fasta.py ' + outReadsFile
    print 'Convert to fasta cmd:' + cmd
    genutils.runCMD(cmd)
    readsFasta = outReadsFile.replace('.fq','.fasta')

    #setup constraint, assume fragments are form 20 to 500 bp
    cmd = 'formcon %s 20 500' % (readsFasta)
    print 'setup constrains: ',cmd
    genutils.runCMDNoFail(cmd)
    
    #cap3 cmd
    cap3Out = readsFasta + '.cap3.screen'
    cmd = 'cap3 %s -i 25 -j 31 -o 16 -s 251 -z 1 -c 10 > %s' % (readsFasta,cap3Out)
    print 'cap3 cmd:',cmd
    genutils.runCMD(cmd)
    
    #get fasta sizes
    contigsFasta = readsFasta + '.cap.contigs'
    fastaSizes = contigsFasta + '.sizes'   

    cmd = 'fastalength %s > %s ' % (contigsFasta,fastaSizes)
    if  os.path.isfile(contigsFasta) is True:
        genutils.runCMD(cmd)


LOG.close()
print 'AND ALL DONE'

