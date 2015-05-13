import sys
import genutils
import brkptalign
import os.path

from optparse import  OptionParser


#############################################################################

USAGE = """python 15_gather-alternate-alleles.py
                                  --in <input of miropeats bp parse >
                                  --out <output bp output table >
                                  --ref <reference genome fasta, must be indexed>
                                  --miropeats_base <miropeats base dir>
                                  --allelesBaseDir <base dir for alleles output>


"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inputFileName', help = 'file of initial miropeats input')
parser.add_option('--out',dest='outputFileName', help = 'file of bp output table')
parser.add_option('--ref',dest='refFasta', help = 'reference genome fasta')
parser.add_option('--miropeats_base',dest='mrDir', help = 'directory for miropeats base')
parser.add_option('--allelesBaseDir',dest='allelesDir', help = 'directory for alleles seq files')


(options,args)=parser.parse_args()

if options.inputFileName is None:
    parser.error('input miropeats result file not given')
if options.outputFileName is None:
    parser.error('output bp result file not given')
if options.refFasta is None:
    parser.error('refFasta file  not given')
if options.mrDir is None:
    parser.error('miropeats dir not given')
if options.allelesDir is None:
    parser.error('allele seqs output dir not given')



    
#############################################################################
if options.mrDir[-1] != '/':
     options.mrDir  += '/'
if options.allelesDir[-1] != '/':
     options.allelesDir  += '/'


refGenomeFasta = options.refFasta
workingDirBase = options.mrDir 

parsedBPFile = options.inputFileName
outResInputfileOut = options.outputFileName

allelesBaseDir = options.allelesDir
bpOutTable = options.outputFileName
parsedBPFile = options.inputFileName



regDelta = 2000
fragmentExtension = 600



inFile = open(parsedBPFile,'r')
outTable = open(bpOutTable,'w')

outTable.write('siteID\tchrom\tInsBP\tgTSDl\talignFragBegin\talignFragEnd\n')
print parsedBPFile
print bpOutTable
print allelesBaseDir

numDid = 0
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    siteID = line[0]
    if siteID == 'siteID':
        continue

    # get all of the info
    data = {}
    brkptalign.populate_all_breakpoint_results(data,line,regDelta,workingDirBase,refGenomeFasta)
 
    print data['siteID']       
    brkptalign.make_alternative_seqs(data,outTable,allelesBaseDir,fragmentExtension)



    numDid += 1
#    if numDid >= 2:
#        break
inFile.close()
outTable.close()
