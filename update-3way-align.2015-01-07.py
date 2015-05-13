import sys
import genutils
import brkptalign
import os.path

refGenomeFasta = '/home/jmkidd/kidd-lab/genomes/canFam3.1/canFam3.1-withUn/canFam3.1.withUn.fa'
assemBaseDir = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/canine-retroseq/SINE/assembly/'
parseOutBaseDir = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/canine-retroseq/SINE/parse_2015-01-06/'
workingDirBase = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/canine-retroseq/SINE/parse_2015-01-06/miropeats/'




outResInputfile = parseOutBaseDir + 'initial-bp.2015-01-06.3wayalign.good.changes'

outResInputfileOut = parseOutBaseDir + 'initial-bp.2015-01-06.3wayalign.good.changes.a2'
#outResInputfileOut = outResInputfile + '.a1'

targetType = 'SINE/tRNA' # for the sliding..

# do we use these things?
regDelta = 2000
minINS = 50
minStart = 30
minStart = 60


inFile = open(outResInputfile,'r')
print outResInputfile
outTable = open(outResInputfileOut,'w')

header = ['siteID','chromosome','ContigName','ContigOreintation','ContigLen','chromFragStart','chromFragEnd','leftChromEnd','leftContigEnd','rightChromStart','rightContigStart']
header.extend(['genomeAlignFragStart','genomeAlignFragEnd','leftFragStart','leftFragEnd','rightFragStart','rightFragEnd'])
header.extend(['leftBPGFrag','leftBPChrom','leftBPContig','rightBPGFrag','rightBPChrom','rightBPContig','TSDlen','insLen'])


header = '\t'.join(header) + '\n'
outTable.write(header)

toCheck = []
#toCheck = ['chr10_22867541','chr10_117703797','chr11_102525453']
#toCheck = ['chr23_52206013']
#toCheck = []
#chromsToDo = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7']
#chromsToDo = ['chr1']
numDid = 0


for line in inFile:
    line = line.rstrip()
#    print '**',line,'**'
    line = line.split('\t')
    siteID = line[0]
    if siteID == 'siteID':
        continue
    if len(toCheck) > 0 and siteID not in toCheck:
        continue
    print line
    print 'YO'
    
    
    data = {}
    data['refGenomeFasta'] = refGenomeFasta
    brkptalign.populate_data_from_mrfile(data,line,regDelta,workingDirBase,addBreaks=True)
    brkptalign.run_rm(data,run=False) #just to get file names, already ran
    brkptalign.populate_previous_align_results(data,line)
    brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = False)  # need to redo becuase of repeated rerurns
    
    
    # now we have some issues here to deal with
    # need to update the coordinates of where we want to be

    if False:
        data['genomeFragAlignBegin'] = int(line[11])
        data['genomeFragAlignEnd'] = int(line[12])
        data['contigLeftFragAlignBegin'] = int(line[13])
        data['contigLeftFragAlignEnd'] = int(line[14])
        data['contigRightFragAlignBegin'] = int(line[15])
        data['contigRightFragAlignEnd'] = int(line[16])



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
    elif changes in ['shift right left','slide right left']: 
        data['contigRightFragAlignBegin'] = data['contigRightFragAlignBegin'] - 3*extendSize
        data['contigRightFragAlignEnd'] = data['contigRightFragAlignEnd'] - 3*extendSize     
        brkptalign.regulate_align_frags(data)
        brkptalign.repeat_3way_align_after_fragment_update(data,updatedBp = True)
    elif changes in ['shift left left','slide left left']: 
        data['contigLeftFragAlignBegin'] = data['contigLeftFragAlignBegin'] - 3*extendSize
        data['contigLeftFragAlignEnd'] = data['contigLeftFragAlignEnd'] - 3*extendSize     
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
            
    # temp manual


    
   # end temp manual

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
    cmd = 'python /home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/canine-retroseq/SINE/scripts/annotate-miropeats-3align.py '
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
#    if numDid > 19:
#        break
#    break
    
inFile.close()
outTable.close()