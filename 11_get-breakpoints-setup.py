import sys
import genutils
import time
import os.path
from optparse import  OptionParser

###############################################################################

USAGE = """11_get-breakpoints-setup.py --in <input file not both gap> --ref <reference genome fasta, must be indexed>
           --assembase <base directory for assembly output>
           --parsebase <base directory for alignment parse>
           

"""
parser = OptionParser(USAGE)
parser.add_option('--ref',dest='refFasta', help = 'reference genome fasta')
parser.add_option('--assembase',dest='assemBaseDir', help = 'base directory for assembly')
parser.add_option('--parsebase',dest='parseBaseDir', help = 'base directory for assembly')
parser.add_option('--in',dest='inFileName', help = 'sel longest fragment')


(options,args)=parser.parse_args()
if options.refFasta is None:
    parser.error('refFasta file  not given')
if options.assemBaseDir is None:
    parser.error('assembly base dir  not given')
if options.parseBaseDir is None:
    parser.error('parse alignment base dir  not given')
if options.inFileName is None:
    parser.error('input file -- not both gap -- not given')

#############################################################################

if options.parseBaseDir[-1] != '/':
    options.parseBaseDir += '/'



###############################################################################
refGenomeFasta = options.refFasta
assemBaseDir = options.assemBaseDir
parseOutBaseDir = options.parseBaseDir

outResFile = options.inFileName + '.getcontig'
print 'Outpint results for contig selection to',outResFile



regDelta = 2000
minINS = 50
minStart = 30
minStart = 60



inFile = open(options.inFileName,'r')
outRES = open(outResFile,'w')
header = ['siteID','contigBest','strand','contigLen','alignLen']
header = '\t'.join(header) + '\n'
outRES.write(header)

numDid = 0
for line in inFile:
    line = line.rstrip()
    line = line.split()
    siteID = line[0]

    print line
    
    contigs = line[2:]
    contigNames = []
    for contig in contigs:
        cName = contig.split(':')[0] 
        contigNames.append(cName)

    # get site ID and setup output dir info
    chrom = siteID.split('_')
    chrom = chrom[0:-1]
    chrom = '_'.join(chrom)

    p = int(siteID.split('_')[-1])
    


    assemDir = assemBaseDir + '/' + chrom + '/' + siteID
    print assemDir

    alignOutDir = parseOutBaseDir + '/' + chrom + '/' + siteID
    cmd = 'mkdir -p ' + alignOutDir
    print cmd
    genutils.runCMD(cmd)
    alignOutDir += '/'

# get genome region
    startBp = p - regDelta
    endBp = p + regDelta
    region = chrom + ':' + str(startBp) + '-' + str(endBp)
    
    genomeFragName = alignOutDir + 'genomeFrag.fa'
    cmd = 'samtools faidx ' + refGenomeFasta + ' ' + region + ' > ' + genomeFragName
    print cmd
    genutils.runCMD(cmd)

    
    
    contigsFileName = assemDir + '/' + 'combined.scaffolds.fa'
    contigSeq = genutils.read_fasta_file_to_list(contigsFileName)
    

    
    hitResults = {}
    for cName in contigNames:
        print 'doing BLAT for',cName
        seq = contigSeq[cName]['seq']        
        contigSeqStr = genutils.add_breaks_to_line(seq)
        contigFileName =   alignOutDir + cName + '.fa'  
        outFile = open(contigFileName,'w')
        outFile.write('>%s\n%s\n' % (cName,contigSeqStr))
        outFile.close()
        
        outPSL = contigFileName + '.blat'        
        cmd = 'blat -tileSize=6 -minScore=15 -oneOff=1 ' + genomeFragName + ' ' + contigFileName + ' ' + outPSL
        print cmd
        genutils.runCMD(cmd)
        
        # parse just the top hit
        hitIn = open(outPSL,'r')
        # check for number of lines in files
        hitLine = ''
        numLines = 0
        hitLines = []
        for line in hitIn:
            numLines += 1
            if numLines >= 6:
                line = line.rstrip()
                line = line.split()
                hitLines.append(line)
        hitIn.close()

#        print 'HITS'
#        print hitLines
        if len(hitLines) > 0:
            maxScore = 0
            maxScore_i = 0
            for i in range(len(hitLines)):
                s = int(hitLines[i][0])
                if s > maxScore:
                    maxScore = s
                    maxScore_i = i
            hitLine = hitLines[maxScore_i]
            print hitLine
            strand = hitLine[8]
            blockSizes = hitLine[18]
            numBlocks = int(hitLine[17])
            if numBlocks == 1:
                hitResults[cName] = {}
                hitResults[cName]['strand'] = strand
                hitResults[cName]['totAlign'] = -1                        
            else:    
                blocks = blockSizes.split(',')
                totBlock = 0
                for i in range(len(blocks)-1):
                    totBlock += int(blocks[i])
                print totBlock
                hitResults[cName] = {}
                hitResults[cName]['strand'] = strand
                hitResults[cName]['totAlign'] = totBlock                        
        else:
            hitResults[cName] = {}
            hitResults[cName]['strand'] = '?'
            hitResults[cName]['totAlign'] = 0                        

        
    # here, did BLAT, so now need to pick a hit
    hitCandidate = []
    m = 0
    
    for cName in contigNames:
        if hitResults[cName]['strand'] == '?':
            continue
        contigLen =  contigSeq[cName]['seqLen']
        totAlign =  hitResults[cName]['totAlign']
        if totAlign == -1:  #flag, means only 1 block hit
            continue
        if (contigLen - totAlign) >= minINS:
            hitCandidate.append([cName,contigLen - totAlign])
            if (contigLen - totAlign) > m:
                m = contigLen - totAlign
    
    numCands = len(hitCandidate)
    print 'have',numCands,'candidates'
    if numCands == 0:
        print 'SKIPPING, have zero candidates'
        nl = [siteID,'NO_CANDIDATE']
        nl = [str(j) for j in nl]
        nl = '\t'.join(nl) + '\n'
        outRES.write(nl)
    else:
        print hitCandidate
        cBest = ''
        for i in hitCandidate:
            if i[1] == m:
                cBest = i[0]
        print 'and best is',cBest,
    
        nl = [siteID,cBest, hitResults[cBest]['strand'],contigSeq[cBest]['seqLen'],hitResults[cBest]['totAlign'] ]
        nl = [str(j) for j in nl]
        nl = '\t'.join(nl) + '\n'
        outRES.write(nl)

    
    
    
    print siteID
    print contigNames
    
    numDid += 1
    
inFile.close()
outRES.close()
