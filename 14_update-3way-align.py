import sys
import genutils
import brkptalign
import os.path

from optparse import  OptionParser


#############################################################################

USAGE = """python 14_update-3way-align.py 
                                  --in <input of miropeats parse with updates>
                                  --out <output for new breakpoints parse>
                                  --ref <reference genome fasta, must be indexed>
                                  --miropeats_base <miropeats base dir>


"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inputFileName', help = 'file of initial miropeats input')
parser.add_option('--out',dest='outputFileName', help = 'file of new miropeats output')
parser.add_option('--ref',dest='refFasta', help = 'reference genome fasta')
parser.add_option('--miropeats_base',dest='mrDir', help = 'directory for temp files')


(options,args)=parser.parse_args()

if options.inputFileName is None:
    parser.error('input miropeats result file not given')
if options.outputFileName is None:
    parser.error('output miropeats result file not given')

if options.refFasta is None:
    parser.error('refFasta file  not given')
if options.mrDir is None:
    parser.error('miropeats dir not given')



    
#############################################################################
if options.mrDir[-1] != '/':
     options.mrDir  += '/'

refGenomeFasta = options.refFasta
workingDirBase = options.mrDir 
outResInputfile = options.inputFileName
outResInputfileOut = options.outputFileName


print 'Writing table of 3 way alignment updated results to',outResInputfileOut


# hardcoded in base scripts dir
# set scriptDir to location where scripts are installed
#scriptDir = 'scripts/'
scriptDir = '../'
regDelta = 2000
targetType = 'SINE/Alu' # for the sliding..


inFile = open(outResInputfile,'r')
outTable = open(outResInputfileOut,'w')

header = ['siteID','chromosome','ContigName','ContigOreintation','ContigLen','chromFragStart','chromFragEnd','leftChromEnd','leftContigEnd','rightChromStart','rightContigStart']
header.extend(['genomeAlignFragStart','genomeAlignFragEnd','leftFragStart','leftFragEnd','rightFragStart','rightFragEnd'])
header.extend(['leftBPGFrag','leftBPChrom','leftBPContig','rightBPGFrag','rightBPChrom','rightBPContig','TSDlen','insLen'])


header = '\t'.join(header) + '\n'
outTable.write(header)

numDid = 0

passed = True

toCheck = []
#toCheck.append('chr1_115543668')
for line in inFile:
    line = line.rstrip()
    line = line.split('\t')
    siteID = line[0]
    if siteID == 'siteID':
        continue
    if len(toCheck) > 0 and siteID not in toCheck:
        continue
    print line
    
    data = {}
    data['refGenomeFasta'] = refGenomeFasta
    brkptalign.populate_data_from_mrfile(data,line,regDelta,workingDirBase,addBreaks=True)
    brkptalign.run_rm(data,run=False) #just to get file names, already ran

    data['genomeFragAlignBegin'] = int(line[11])
    data['genomeFragAlignEnd'] = int(line[12])
    data['contigLeftFragAlignBegin'] = int(line[13])
    data['contigLeftFragAlignEnd'] = int(line[14])
    data['contigRightFragAlignBegin'] = int(line[15])
    data['contigRightFragAlignEnd'] = int(line[16])


    brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True) # because things did not match
    brkptalign.populate_previous_align_results(data,line)

    print siteID
    extendSize = 25

    if len(line) == 29:
        print 'UPDATES'
        leftContigMP = line[25]
        rightContigMP = line[26]
        genomeMP = line[27]
        changes = line[28]
    else:
        changes = 'NO_CHANGE'
    
    if changes == 'double':
        changes = 'NO_CHANGE'
    
    changes = changes.lower() # fix them to have less
    if changes in ['change','changes']:
        print 'manual contig extraction updates'
        if leftContigMP != '':
            leftContigMP = int(leftContigMP)
            print 'leftMP',leftContigMP
            data['contigLeftFragAlignBegin'] = leftContigMP - 3*extendSize
            data['contigLeftFragAlignEnd'] = leftContigMP + 3*extendSize
        if rightContigMP != '':
            rightContigMP = int(rightContigMP)
            print 'rightMP',rightContigMP
            print data['contigRightFragAlignBegin'],data['contigRightFragAlignEnd']
            data['contigRightFragAlignBegin'] = rightContigMP - 3*extendSize
            data['contigRightFragAlignEnd'] = rightContigMP + 3*extendSize
        if genomeMP != '':
            genomeMP = int(genomeMP)
            print 'genomeMP',genomeMP
            data['genomeFragAlignBegin'] = data['chromFragStart'] - 1 + genomeMP - 6*extendSize
            data['genomeFragAlignEnd'] = data['chromFragStart'] - 1 + genomeMP + 6*extendSize            
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['center chrom','recenter chrom']:        
        print 'center chrom to',data['leftBpChromCoords']
        data['genomeFragAlignBegin'] = data['leftBpChromCoords']- 6*extendSize
        data['genomeFragAlignEnd'] = data['leftBpChromCoords']  + 6*extendSize            
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['NO_CHANGE','no_change']:
        print 'NO CHANGE, still do update'
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)
    elif changes in ['extend-all','extend all', 'expand all']: 
        data['contigRightFragAlignBegin'] = data['contigRightFragAlignBegin'] - 2*extendSize
        data['contigRightFragAlignEnd'] = data['contigRightFragAlignEnd'] + 2*extendSize     
        data['contigLeftFragAlignBegin'] = data['contigLeftFragAlignBegin'] - 2*extendSize
        data['contigLeftFragAlignEnd'] = data['contigLeftFragAlignEnd'] + 2*extendSize     
        data['genomeFragAlignBegin'] = data['genomeFragAlignBegin'] - 2*extendSize
        data['genomeFragAlignEnd'] = data['genomeFragAlignEnd'] + 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)
    elif changes in ['extend right']: 
        data['contigRightFragAlignBegin'] = data['contigRightFragAlignBegin'] - 2*extendSize
        data['contigRightFragAlignEnd'] = data['contigRightFragAlignEnd'] + 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)
    elif changes in ['extend chrom','extend genome','enlarge chromosome','expand genome','expand chrom']:        
        print 'Extend Chrom!'
        data['genomeFragAlignBegin'] = data['genomeFragAlignBegin'] - 2*extendSize
        data['genomeFragAlignEnd'] = data['genomeFragAlignEnd'] + 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['expand genome a lot']:        
        print 'Extend Chrom!'
        data['genomeFragAlignBegin'] = data['genomeFragAlignBegin'] - 6*extendSize
        data['genomeFragAlignEnd'] = data['genomeFragAlignEnd'] + 6*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['extend genome right']:        
        print 'Extend Chrom!'
        data['genomeFragAlignEnd'] = data['genomeFragAlignEnd'] + 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['extend genome left','expand genome left']:        
        print 'Extend Chrom!'
        data['genomeFragAlignBegin'] = data['genomeFragAlignBegin'] - 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['shorten left']:        
        print 'shorten left'
        data['contigLeftFragAlignEnd'] = data['contigLeftFragAlignEnd'] - 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['shorten right']:        
        print 'shorten right'
        data['contigRightFragAlignBegin'] = data['contigRightFragAlignBegin'] + 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['extend right right']:        
        print 'extend right right'
        data['contigRightFragAlignEnd'] = data['contigRightFragAlignEnd'] + 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['extend right left']:        
        print 'extend right left'
        data['contigRightFragAlignBegin'] = data['contigRightFragAlignBegin'] - 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['slide left left and slide genome right']:        
        print 'extend right right'
        data['genomeFragAlignBegin'] = data['genomeFragAlignBegin'] + 2*extendSize
        data['genomeFragAlignEnd'] = data['genomeFragAlignEnd'] + 2*extendSize     
        data['contigLeftFragAlignBegin'] = data['contigLeftFragAlignBegin'] - 2*extendSize
        data['contigLeftFragAlignEnd'] = data['contigLeftFragAlignEnd'] - 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['slide genome left']:        
        print 'extend right right'
        data['genomeFragAlignBegin'] = data['genomeFragAlignBegin'] - 2*extendSize
        data['genomeFragAlignEnd'] = data['genomeFragAlignEnd'] - 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['extend left left']:        
        print 'extend left left'
        data['contigLeftFragAlignBegin'] = data['contigLeftFragAlignBegin'] - 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['extend left right']:        
        print 'extend left left'
        data['contigLeftFragAlignEnd'] = data['contigLeftFragAlignEnd'] + 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['extend left left and right right']:        
        print 'extend left left and right right'
        data['contigLeftFragAlignBegin'] = data['contigLeftFragAlignBegin'] - 2*extendSize
        data['contigRightFragAlignEnd'] = data['contigRightFragAlignEnd'] + 2*extendSize     
             
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)        
    elif changes in ['slide right right']: 
        data['contigRightFragAlignBegin'] = data['contigRightFragAlignBegin'] + 2*extendSize
        data['contigRightFragAlignEnd'] = data['contigRightFragAlignEnd'] + 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)
    elif changes in ['slide genome right']: 
        data['genomeFragAlignBegin'] = data['genomeFragAlignBegin'] + 2*extendSize
        data['genomeFragAlignEnd'] = data['genomeFragAlignEnd'] + 2*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)
    elif changes in ['extend-genome-check']: 
        print 'genome check and extend'
        if (data['genomeFragAlignEnd'] - data['genomeFragAlignBegin'] ) >= 1000:
            print 'genomeMP',genomeMP
            data['genomeFragAlignBegin'] = data['chromFragStart'] - 1 + 2000 - 6*extendSize
            data['genomeFragAlignEnd'] = data['chromFragStart'] - 1 + 2000 + 6*extendSize            
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)
    else:
        print '??'
        print data['siteID'],changes
        print '??'
        sys.exit()

    # ok, ready now to update the rest
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
inFile.close()
outTable.close()