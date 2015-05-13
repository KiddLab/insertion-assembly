import os
import sys
import genutils
from optparse import  OptionParser




##############################################################################################################################################################
def process_results(regionName,outDir,outFile):
        
    contigsFile = outDir + '/combined.scaffolds.fa'
    if os.path.isfile(contigsFile) is True:
        fastaFile = contigsFile
        rmFile = fastaFile + '.out'
    else:
        outFile.write('%s\tNO_CONTIGS\t%i\n' % (regionName,0 ))
        return            
    if os.path.isfile(rmFile) is False:
        outFile.write('%s\tNO_RMOUT\t%i\n' % (regionName,0 ))
        return
            
    # get list of fastaFiles
    fastaDict = genutils.read_fasta_file_to_list(fastaFile)
    rmLines = genutils.read_rm_file(rmFile)
    seqToPrint = {}
    hits = []
    for rm in rmLines:
        # change this to match what you are looking for!
        if rm[9][:3] == 'Alu':
            seqName = rm[4]
            start = int(rm[5])
            end = int(rm[6])
            before = start -1
            after = fastaDict[seqName]['seqLen'] - end

            report = '%s:%i:%i:%i-%i' % (seqName,before,after,start,end) 
            hits.append(report)
            seqToPrint[seqName] = 1
    hJ = '\t'.join(hits)
    outFile.write('%s\tHAS_SEQ\t%i\t%s\n' % (regionName,len(hits),hJ))   
###############################################################################
#############################################################################

USAGE = """06_parse-sine-cap3-only-assembly.py --vcf <input vcf file> --outdirbase <base directory for assembly output>
                                               --out <assembly summary file name>

"""
parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcfFile', help = 'input vcf file name')
parser.add_option('--outdirbase',dest='outDirBase', help = 'base directory for assembly output')
parser.add_option('--out',dest='outFileName', help = 'assembly output summary file name')


(options,args)=parser.parse_args()

if options.vcfFile is None:
    parser.error('vcf file name not given')
if options.outDirBase is None:
    parser.error('base dir for assembly output not given')
if options.outFileName is None:
    parser.error('output file name for parsed assembly table not given')
    
#############################################################################


outDirBase = options.outDirBase
inputVCFCalls = options.vcfFile

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
    regionsToDo.append([c,k])
inFile.close()

print 'Read in %i regions to do' % (len(regionsToDo))

outFile = open(options.outFileName,'w')

n = 0
for regionName in regionsToDo:
    outDir = outDirBase + '/' + regionName[0] + '/' + regionName[1]
    
    if os.path.isdir(outDir) is False:
        outFile.write('%s\tNO_CAP3_directory\n' % (regionName[1]))
    
    process_results(regionName[1],outDir,outFile)
    n += 1
    if n%1000==0:
        print 'Doing %i of %i ...' % (n,len(regionsToDo))
        
print 'ALL DONE'
    

outFile.close()


