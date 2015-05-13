import genutils
import os
import sys
from optparse import  OptionParser


#############################################################################

USAGE = """combine-contigs-based-on-links.py --vcf <input vcf file> --outdirbase <base directory for assembly output>

        Combines contigs into scaffolds (based on links) and runs RepeatMasker
  
"""
parser = OptionParser(USAGE)
parser.add_option('--vcf',dest='vcfFile', help = 'input vcf file name')
parser.add_option('--outdirbase',dest='outDirBase', help = 'base directory for assembly output')

(options,args)=parser.parse_args()

if options.vcfFile is None:
    parser.error('vcf file name not given')
if options.outDirBase is None:
    parser.error('base dir for assembly output not given')

    
#############################################################################
inputVCFCalls = options.vcfFile
outDirBase = options.outDirBase 


intervalsToReadNames = {}
nameToSeq = {}

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
    regionsToDo.append([k,c])
inFile.close()

print 'Read in %i regions to do' % (len(regionsToDo))

gapSeq = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'

for region in regionsToDo:
    c = region[1]
    region = region[0]

    print '***',region,'***'
    outPutDir = outDirBase + '/' + c + '/' + region 
    outPutDir += '/'
    
    print region
    print outPutDir
    
    if os.path.isdir(outPutDir) is False:
        print 'Dir not exist, skipping'
        continue
    
    contigsFile = outPutDir + 'reads.fasta.cap.contigs'
    contigFileSize = os.path.getsize(contigsFile)
    if contigFileSize == 0:
        newContigFile = outPutDir + 'combined.scaffolds.fa'
        outFile = open(newContigFile,'w')
        outFile.close()
        continue
        
    contigs = genutils.read_fasta_file_to_list(contigsFile)
    
    linksFile = outPutDir + 'reads.fasta.cap.contigs.links'
    conLinks = []
    inFile = open(linksFile,'r')
    for line in inFile:
        line = line.rstrip()
        if line == '':
            continue
        line = line.split(' ')
        c1 = line[0]
        c2 = line[2]    
        conLinks.append([c1,c2])
    inFile.close()
    conLinks = conLinks[::-1]
    print 'have %i contig links' % len(conLinks)
    for i in conLinks:
        contigs[i[0]]['seq'] += gapSeq + contigs[i[1]]['seq']
        del contigs[i[1]]
    
    newContigFile = outPutDir + 'combined.scaffolds.fa'
    outFile = open(newContigFile,'w')
    cNames = contigs.keys()
    cNames.sort()
    for c in cNames:
        s = genutils.add_breaks_to_line(contigs[c]['seq'])
        cN = c.replace('Contig','Scafftig')
        outFile.write('>%s\n%s\n' % (cN,s))    
    outFile.close()

    fastaSizes = newContigFile + '.sizes'   

    cmd = 'fastalength %s > %s ' % (newContigFile,fastaSizes)
    genutils.runCMD(cmd)
    # change this based on species
    cmd = 'RepeatMasker --species human ' + newContigFile
#    cmd = 'RepeatMasker --species dog ' + newContigFile    
    print cmd
    genutils.runCMD(cmd)


print 'ALL DONE'
