import sys
import genutils
import brkptalign
import os.path

# starts with trying to get breakpoint from miropeats
from optparse import  OptionParser
###############################################################################

USAGE = """12_get-breakpoints-initialMR.py  --ref <reference genome fasta, must be indexed> --tmp <tmp dir>
               --assembase <base directory for assembly output>
               --miropeats_base <miropeats base dir>
               --in <input file with getcontig> 
              --out <output file for initial miropeats results>
           
"""
parser = OptionParser(USAGE)
parser.add_option('--ref',dest='refFasta', help = 'reference genome fasta')
parser.add_option('--tmp',dest='tmpDir', help = 'directory for temp files')
parser.add_option('--assembase',dest='assemBaseDir', help = 'base directory for assembly')
parser.add_option('--miropeats_base',dest='mrDir', help = 'directory for temp files')
parser.add_option('--in',dest='inFileName', help = 'input file with getcontig results')
parser.add_option('--out',dest='outFileName', help = 'out file for initial miropeats results')



(options,args)=parser.parse_args()
if options.refFasta is None:
    parser.error('refFasta file  not given')
if options.tmpDir is None:
    parser.error('tmp dir not given')
if options.mrDir is None:
    parser.error('mrDir dir not given')
if options.assemBaseDir is None:
    parser.error('assembly base dir  not given')
if options.inFileName is None:
    parser.error('input file name with get contig results  not given')
if options.outFileName is None:
    parser.error('output file name with initial miropeats results  not given')
#############################################################################


if options.assemBaseDir[-1] != '/':
     options.assemBaseDir  += '/'

if options.mrDir[-1] != '/':
     options.mrDir  += '/'

if options.tmpDir[-1] != '/':
     options.tmpDir  += '/'


refGenomeFasta = options.refFasta
workingDirBase = options.mrDir
assemBaseDir =  options.assemBaseDir
outResFile = options.inFileName
outResFile2 = options.outFileName

# hardcoded in base scripts dir
# set scriptDir to location where scripts are installed
#scriptDir = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/RetroSeq-HGDP/retroseq-alu-additional/results/assembly/hg18/scripts/'
scriptDir = '../'


regDelta = 2000
minINS = 50
minStart = 30
minStart = 60


inFile = open(outResFile,'r')
outTable = open(outResFile2,'w')

header = ['siteID','chromosome','ContigName','ContigOreintation','ContigLen','chromFragStart','chromFragEnd','leftChromEnd','leftContigEnd','rightChromStart','rightContigStart','chromTSDSize','contigInsSize']
header = '\t'.join(header) + '\n'
outTable.write(header)


numDid = 0
for line in inFile:
    line = line.rstrip()
    line = line.split()
    siteID = line[0]
    if siteID == 'siteID':
        continue
    if line[1] == 'NO_CANDIDATE':
        continue

    print line
    
    data = {}
    data['tmpDir'] = options.tmpDir
    if os.path.isdir(data['tmpDir']) is False:
        cmd = 'mkdir ' + data['tmpDir']
        print 'making tmp dir'
        print cmd
        genutils.runCMD(cmd)
    
    data['refGenomeFasta'] = refGenomeFasta
    
    data['siteID'] = siteID

    chrom = siteID.split('_')
    chrom = chrom[0:-1]
    chrom = '_'.join(chrom)

    p = int(siteID.split('_')[-1])
    startBp = p - regDelta
    endBp = p + regDelta
    data['chromName'] = chrom
    data['chromFragStart'] =startBp
    data['chromFragEnd'] = endBp
    
    data['contigName'] = line[1]
    data['contigDir'] = line[2]


    originalContigsDir = assemBaseDir + data['chromName'] + '/' + data['siteID'] 
    data['originalContigsDir']  = originalContigsDir


    alignOutDir = workingDirBase + data['chromName'] + '/' + data['siteID'] 
    cmd = 'mkdir -p ' + alignOutDir
    genutils.runCMD(cmd)
    data['alignOutDir'] = alignOutDir


    brkptalign.get_contig_seq(data)
    brkptalign.get_genome_frag(data)
    
    brkptalign.get_contig_gaps(data,run=True)
    

    nl = [data['siteID'],data['chromName'],data['contigName'],data['contigDir']]
    nl.append(str(data['contigLen']))

    nl.append(str(data['chromFragStart']))
    nl.append(str(data['chromFragEnd']))

    brkptalign.run_miropeats(data)
    brkptalign.parse_miropeats_hits(data)
    nl.append(str(data['leftChromFragEnd']))
    nl.append(str(data['leftContigEnd']))

    nl.append(str(data['rightChromFragStart']))
    nl.append(str(data['rightContigStart']))

    if data['leftChromFragEnd'] == 'NA':
        data['chromTSDSize'] = 'NA'
        data['contigInsSize'] = 'NA'    
    else:
        data['chromTSDSize'] = data['leftChromFragEnd'] - data['rightChromFragStart'] + 1
        data['contigInsSize'] = data['rightContigStart'] - data['leftContigEnd'] + 1
        
    nl.append(str(data['chromTSDSize']))
    nl.append(str(data['contigInsSize']))

    brkptalign.run_rm(data,run=True)

    # run initial parse make image
    
    if data['miroOutPS'] != 'FAILURE':    
		cmd = 'python ' + scriptDir + 'annotate-miropeats.py '    
		cmd += ' --miroin ' + data['miroOutPS']
		cmd += ' --topRM ' + data['genomeFragFileNameRM']
		cmd += ' --bottomRM ' + data['contigSeqFileNameRM']
		cmd += ' --clonebreak ' + str(data['rightContigStart']) + ',' + str(data['leftContigEnd'])
		cmd += ' --chrombreak ' + str(data['leftChromFragEnd']) + ',' + str(data['rightChromFragStart'])
		cmd += ' --siteID ' + data['siteID']
		cmd += ' --tsd ' + str(data['chromTSDSize'])
		cmd += ' --ins ' + str(data['contigInsSize'])
		cmd += ' --clonegap ' + data['contigSeqGapsFileName']
	
		print cmd
		genutils.runCMD(cmd)

		data['psAnnotated'] = data['miroOutPS'] + '.annotated.ps'
		data['psAnnotatedPDF'] = data['alignOutDir'] + '/' + data['siteID'] + '.initialMR.annotated.pdf'
	
		if os.path.isfile(data['psAnnotatedPDF']) is True:
			cmd = 'rm ' + data['psAnnotatedPDF']
			genutils.runCMD(cmd)
	
		cmd = 'ps2pdf %s %s' % (data['psAnnotated'],data['psAnnotatedPDF'])
		print cmd
		genutils.runCMD(cmd)
    


    nl = '\t'.join(nl) + '\n'
    outTable.write(nl)


inFile.close()
outTable.close()