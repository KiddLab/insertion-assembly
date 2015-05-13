import sys
import genutils
import brkptalign
import os.path
from optparse import  OptionParser


#############################################################################

USAGE = """python 13_get-breakpoints-3way-align.py  
                                  --in <input of initial miropeats parse>
                                  --ref <reference genome fasta, must be indexed>
                                  --miropeats_base <miropeats base dir>


"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inputFileName', help = 'file of initial miropeats input')

parser.add_option('--ref',dest='refFasta', help = 'reference genome fasta')
parser.add_option('--miropeats_base',dest='mrDir', help = 'directory for temp files')


(options,args)=parser.parse_args()

if options.inputFileName is None:
    parser.error('input miropeats result file not given')
if options.refFasta is None:
    parser.error('refFasta file  not given')
if options.mrDir is None:
    parser.error('miropeats dir not given')



    
#############################################################################
if options.mrDir[-1] != '/':
     options.mrDir  += '/'

refGenomeFasta = options.refFasta
workingDirBase = options.mrDir 

# hardcoded in base scripts dir
# set scriptDir to location where scripts are installed
#scriptDir = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/RetroSeq-HGDP/retroseq-alu-additional/results/assembly/hg18/scripts/'
scriptDir = '../'
regDelta = 2000


outResFile = options.inputFileName
outResFile2 = outResFile + '.3wayalign'
print 'Writing table of 3 way alignment results to',outResFile2

#assemBaseDir = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/canine-retroseq/SINE/assembly/'
#parseOutBaseDir = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/canine-retroseq/SINE/parse_2015-01-06/'
#workingDirBase = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/canine-retroseq/SINE/parse_2015-01-06/miropeats/'




inFile = open(outResFile,'r')
outTable = open(outResFile2,'w')

header = ['siteID','chromosome','ContigName','ContigOreintation','ContigLen','chromFragStart','chromFragEnd','leftChromEnd','leftContigEnd','rightChromStart','rightContigStart']
header.extend(['genomeAlignFragStart','genomeAlignFragEnd','leftFragStart','leftFragEnd','rightFragStart','rightFragEnd'])
header.extend(['leftBPGFrag','leftBPChrom','leftBPContig','rightBPGFrag','rightBPChrom','rightBPContig','TSDlen','insLen'])


header = '\t'.join(header) + '\n'
outTable.write(header)

numDid = 0
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    siteID = line[0]
    if siteID == 'siteID':
        continue    
    print line
    #setup dictionary of data for each site
    data = {}
    data['refGenomeFasta'] = refGenomeFasta
    brkptalign.populate_data_from_mrfile(data,line,regDelta,workingDirBase,addBreaks=True)
    brkptalign.run_rm(data,run=False) #just to get file names, already ran

    brkptalign.run_3way_align(data)


    # setup output to print
    nl = [data['siteID'],data['chromName'],data['contigName'],data['contigDir']]
    nl.append(str(data['contigLen']))
    nl.append(str(data['chromFragStart']))
    nl.append(str(data['chromFragEnd']))
    nl.append(str(data['leftChromFragEnd']))
    nl.append(str(data['leftContigEnd']))
    nl.append(str(data['rightChromFragStart']))
    nl.append(str(data['rightContigStart']))


    nl.append(data['genomeFragAlignBegin'])
    nl.append(data['genomeFragAlignEnd'])

    nl.append(data['contigLeftFragAlignBegin'])
    nl.append(data['contigLeftFragAlignEnd'])

    nl.append(data['contigRightFragAlignBegin'])
    nl.append(data['contigRightFragAlignEnd'])


    nl.append(data['leftBpGenomeFragCoords'])
    nl.append(data['leftBpChromCoords'])
    nl.append(data['leftBpContigCoord'])

    nl.append(data['rightBpGenomeFragCoords'])
    nl.append(data['rightBpChromCoords'])
    nl.append(data['rightBpContigCoord'])


    gTSDs = data['rightBpChromCoords']
    gTSDe = data['leftBpChromCoords']
    gTSDl = gTSDe - gTSDs + 1
    nl.append(str(gTSDl))

    insStart = data['leftBpContigCoord'] + 1
    insEnd = data['rightBpContigCoord'] - 1

    s = insEnd-insStart+1
    nl.append(str(s))
    
    data['gTSDl'] = gTSDl
    data['insLen'] = s
    
    dataPickleFileName = data['alignOutDir'] + '/' + data['siteID'] + '.data.pickle'
    data['dataPickleFileName'] = dataPickleFileName
    brkptalign.write_pickle_dictionary(data,dataPickleFileName)

#    brkptalign.print_dictionary_keys(data)

    
    nl = [str(j) for j in nl]
    nl = '\t'.join(nl) + '\n'
    outTable.write(nl)
    
    # draw the combined align'
    cmd = 'python '  + scriptDir + 'annotate-miropeats-3align.py '
    cmd += ' --pickle ' + dataPickleFileName
    print cmd
    genutils.runCMD(cmd)

    data['annotatedBPPS'] = data['miroOutPS'] + '.breakpoints.annotated.ps'
    data['psBPAnnotatedPDF'] = data['alignOutDir'] + '/' + data['siteID'] + '.3way.combined.annotated.pdf'

    if os.path.isfile(data['psBPAnnotatedPDF']) is True:
        cmd = 'rm ' + data['psBPAnnotatedPDF']
        genutils.runCMD(cmd)

    cmd = 'ps2pdf %s %s' % (data['annotatedBPPS'],data['psBPAnnotatedPDF'])
    print cmd
    genutils.runCMD(cmd)

    numDid += 1
#    if numDid > 1:
#        break
    

inFile.close()
outTable.close()