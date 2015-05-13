# brkptalign.py
# a set of routines that we want to use for parsing and processing breakpoints and alignments

import genutils
import sys
import os.path
import glob
import pickle

###############################################################################
def write_pickle_dictionary(data,fileName):
    myFile = open(fileName,'w')
    pickle.dump(data,myFile)
    myFile.close()
###############################################################################
def read_pickle_dictionary(fileName):
    myFile = open(fileName,'r')
    myData = pickle.load(myFile)
    myFile.close()
    return myData
###############################################################################
def print_dictionary(d):
    k = d.keys()
    k.sort()
    for i in k:
        print i,d[i]
###############################################################################
def print_dictionary_keys(d):
    k = d.keys()
    k.sort()
    for i in k:
        print i
###############################################################################
def get_contig_seq(data):
    print_dictionary(data)
    contigsFileName = data['originalContigsDir'] + '/' + 'combined.scaffolds.fa'
    contigSeq = genutils.read_fasta_file_to_list(contigsFileName)
    
    cSeq = contigSeq[data['contigName']]['seq']
    if data['contigDir'] == '-':
        cSeq = genutils.revcomp(cSeq)
        
    data['contigSeqFileName'] = data['alignOutDir'] + '/' + 'Contig.genomedir.fa'
    data['contigSeqGenomeDir'] = cSeq
    contigSeqStr = genutils.add_breaks_to_line(cSeq)    
    outFile = open(data['contigSeqFileName'],'w')
    outFile.write('>contig\n%s\n' % contigSeqStr)
    outFile.close()
    
    contigSeq = genutils.read_fasta_to_string(data['contigSeqFileName'])
    data['contigLen'] = len(data['contigSeqGenomeDir'])       
###############################################################################
def get_genome_frag(data):
    region = data['chromName'] + ':' + str(data['chromFragStart']) + '-' + str(data['chromFragEnd'])
    
    data['genomeFragFileName'] = data['alignOutDir'] + '/' + 'genomeFrag.fa'
    
    cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + data['genomeFragFileName']
    print cmd
    genutils.runCMD(cmd)
    
    genomeSeq = genutils.read_fasta_to_string(data['genomeFragFileName'])
    genomeSeq = genomeSeq.upper()
    data['genomeFragSeq'] = genomeSeq
    
###############################################################################
def run_miropeats(data):
    if 'miropeatSValue' in data:
        s = data['miropeatSValue']
    else:
        s = 80
        s = 40 # for the dogs...
        data['miropeatSValue'] = s
    
    
        
    data['miroOutPS'] = data['alignOutDir'] + '/' + 'miropeats.' + str(s) + '.ps'
    data['miroOutInfo'] = data['alignOutDir'] + '/' + 'miropeats.' + str(s) + '.out'
    
    if 'tmpDir' in data:
        tmpDir = data['tmpDir']
    else:
        tmpDir = '/home/jmkidd/kidd-lab-scratch/jmkidd-projects/tmp/'
    
    tmpGenome = tmpDir + 'genome.fa'
    tmpContig = tmpDir + 'contig.fa'
    tmpMRPS = tmpDir + 'tmp.MRPS'
    tmpMROUT = tmpDir + 'tmp.MROUT'
    
    cmd = 'cp %s %s' % (data['genomeFragFileName'],tmpGenome)
    print cmd
    genutils.runCMD(cmd)

    cmd = 'cp %s %s' % (data['contigSeqFileName'],tmpContig)
    print cmd
    genutils.runCMD(cmd)
    
    cmd = 'miropeats -s %i -onlyinter -o %s  -seq %s  -seq %s > %s'  % (s,tmpMRPS,tmpGenome,tmpContig,tmpMROUT)
    print cmd
    genutils.runCMD(cmd)
    
    #cp
    if os.path.isfile(tmpMRPS) is True:    
		cmd = 'cp %s %s' % (tmpMRPS,data['miroOutPS'])
		print cmd
		genutils.runCMD(cmd)
    else:
	    data['miroOutPS'] = 'FAILURE'

    cmd = 'cp %s %s' % (tmpMROUT,data['miroOutInfo'])
    print cmd
    genutils.runCMD(cmd)
    
    # clean up
    
    cmd = 'rm %s %s' % (tmpGenome,tmpContig)
    print cmd
    genutils.runCMD(cmd)
    
    if os.path.isfile(tmpMRPS) is True:
         cmd = 'rm ' + tmpMRPS
         print cmd
         genutils.runCMD(cmd)

    if os.path.isfile(tmpMROUT) is True:
         cmd = 'rm ' + tmpMROUT
         print cmd
         genutils.runCMD(cmd)

  
###############################################################################
# try to guess hit from miropeats..
def parse_miropeats_hits(data):
    print data['miroOutInfo']
    inFile = open(data['miroOutInfo'],'r')
    atHits = False
    hitLines = []
    while True:
        line = inFile.readline()
        if line == '':
            break
        if line[0:9] == '## Sorted':
            atHits = True
            continue
        if atHits is True:
            line = line.rstrip()
            line = line.split()
            # assume is contig first
            l = [int(line[2]),int(line[3]),int(line[7]),int(line[8])]
            hitLines.append(l)
    inFile.close()

    hitLines = merge_mr_hits(hitLines)

    if len(hitLines) == 0:
		data['leftChromFragEnd'] = 'NA'
		data['leftContigEnd'] = 'NA'
		data['rightChromFragStart'] = 'NA'
		data['rightContigStart'] = 'NA'
    elif len(hitLines) == 1:
        fosmidLeftUnMatched = hitLines[0][0] - 1
        fosmidRightUnMatched = data['contigLen'] - hitLines[0][1]
        if fosmidLeftUnMatched > fosmidRightUnMatched: #unmatched left
            data['rightChromFragStart'] = hitLines[0][2]
            data['rightContigStart'] = hitLines[0][0]
            data['leftChromFragEnd'] = data['rightChromFragStart'] - 1
            data['leftContigEnd'] = 1
        else: # umatched right
            data['leftChromFragEnd'] = hitLines[0][3]
            data['leftContigEnd'] = hitLines[0][1]
            data['rightChromFragStart'] = data['leftChromFragEnd'] + 1
            data['rightContigStart'] = data['contigLen']
    else:  # > 1
        contigStarts = []
        contigEnds = []
        for i in hitLines:
            contigStarts.append(i[0])
            contigEnds.append(i[1])
        contigStarts.sort()
        contigEnds.sort()        
        # go through and get end of first one from left
        # we are doing a very simplistic parsing routine...
        for i in hitLines:
            if int(i[0]) == contigStarts[0]:
                data['leftChromFragEnd'] = i[3]
                data['leftContigEnd'] = i[1]
            if int(i[1]) == contigEnds[-1]:
                data['rightChromFragStart'] = i[2]
                data['rightContigStart'] = i[0]
###############################################################################
def parse_miropeats_hits_useRM(data):
    targetType = 'SINEC_Cf'
    print data['miroOutInfo']
    inFile = open(data['miroOutInfo'],'r')
    atHits = False
    hitLines = []
    while True:
        line = inFile.readline()
        if line == '':
            break
        if line[0:9] == '## Sorted':
            atHits = True
            continue
        if atHits is True:
            line = line.rstrip()
            line = line.split()
            # assume is contig first
            l = [int(line[2]),int(line[3]),int(line[7]),int(line[8])]
            hitLines.append(l)
    inFile.close()
    print hitLines

    # might not want to do this
    hitLines = merge_mr_hits(hitLines)
    print hitLines
    
    
    rmLines = read_rm_file(data['contigSeqFileNameRM'])
    rmStarts = []
    rmEnds = []
    for i in rmLines:    
        if i[9][:8] != targetType:
            continue
        rStart = int(i[5])
        rmStarts.append(rStart)
        rEnd = int(i[6])
        rmEnds.append(rEnd)
    rmStarts.sort()
    rmEnds.sort()

    if len(hitLines) == 0:
		data['leftChromFragEnd'] = 'NA'
		data['leftContigEnd'] = 'NA'
		data['rightChromFragStart'] = 'NA'
		data['rightContigStart'] = 'NA'
    elif len(hitLines) == 1:
        fosmidLeftUnMatched = hitLines[0][0] - 1
        fosmidRightUnMatched = data['contigLen'] - hitLines[0][1]
        if fosmidLeftUnMatched > fosmidRightUnMatched: #unmatched left
            data['rightChromFragStart'] = hitLines[0][2]
            data['rightContigStart'] = hitLines[0][0]
            data['leftChromFragEnd'] = data['rightChromFragStart'] - 1
            data['leftContigEnd'] = 1
        else: # umatched right
            data['leftChromFragEnd'] = hitLines[0][3]
            data['leftContigEnd'] = hitLines[0][1]
            data['rightChromFragStart'] = data['leftChromFragEnd'] + 1
            data['rightContigStart'] = data['contigLen']
    else:  # > 1
        print 'RMs'
        print rmStarts
        print rmEnds
        endTarget = hitLines[0][0]

        startTarget = hitLines[0][1]

        startDelta = 100000
        endDelta = 100000
        for i in hitLines:
             if abs(i[0] - rmEnds[-1]) < endDelta:
                 endTarget = i[0]
                 endDelta = abs(endTarget - rmEnds[-1])

             if abs(i[1] - rmStarts[0]) < startDelta:
                 startTarget = i[1]
                 startDelta = abs(startTarget - rmStarts[0])
                
        print 'startTarget',startTarget
        print 'endTarget',endTarget

        # go through and get one that flanks the insertion        
        for i in hitLines:
            if i[0] == endTarget:
                print 'match end'
                data['rightChromFragStart'] = i[2]
                data['rightContigStart'] = i[0]
            if i[1] == startTarget:
                print 'match start'
                data['leftChromFragEnd'] = i[3]
                data['leftContigEnd'] = i[1]
###############################################################################

def run_rm(data,run=True):
    # change species here
    cmd = 'RepeatMasker ' + data['genomeFragFileName']

    if run is True:
        genutils.runCMD(cmd)

    data['genomeFragFileNameRM'] = data['genomeFragFileName'] + '.out'

    # change species here
    cmd = 'RepeatMasker  ' + data['contigSeqFileName']
#    print cmd
    if run is True:
        genutils.runCMD(cmd)

    data['contigSeqFileNameRM'] = data['contigSeqFileName'] + '.out'    
###############################################################################
def get_contig_gaps(data,run=True):
    data['contigSeqGapsFileName'] = data['contigSeqFileName'] + '.gaps'
    cmd = 'get_gaps.pl ' + data['contigSeqFileName'] + ' > ' + data['contigSeqGapsFileName']
    if run is True:
        genutils.runCMD(cmd)    
    get_genome_gaps(data,run)
###############################################################################
def get_genome_gaps(data,run=True):
    data['genomeFragGapsFileName'] = data['genomeFragFileName'] + '.gaps'
    cmd = 'get_gaps.pl ' + data['genomeFragFileName'] + ' > ' + data['genomeFragGapsFileName']
    if run is True:
        genutils.runCMD(cmd)    
###############################################################################
# populate data dictionary from file of initial mr output
def populate_data_from_mrfile(data,line,regDelta,workingDirBase,addBreaks=False):
    data['siteID'] = line[0]
    data['chromName'] = line[1]
    data['contigName'] = line[2]
    data['contigDir'] = line[3]

    data['chromFragStart'] = int(line[5])
    data['chromFragEnd'] = int(line[6])    

    alignOutDir = workingDirBase + data['chromName'] + '/' + data['siteID'] 
    data['alignOutDir'] = alignOutDir


    # contig seq
    data['contigSeqFileName'] = data['alignOutDir'] + '/' + 'Contig.genomedir.fa'
    contigSeq = genutils.read_fasta_to_string(data['contigSeqFileName'])
    data['contigSeqGenomeDir'] = contigSeq
    data['contigLen'] = len(data['contigSeqGenomeDir'])
    if data['contigLen'] != int(line[4]):
        print 'Contig len mismatch!'
        print line
        print data['contigLen']
        sys.exit()
    # genome seq       
    data['genomeFragFileName'] = data['alignOutDir'] + '/' + 'genomeFrag.fa'    
    genomeSeq = genutils.read_fasta_to_string(data['genomeFragFileName'])
    genomeSeq = genomeSeq.upper()
    data['genomeFragSeq'] = genomeSeq

    # miropeats
#    s = 80
# decided to use 40 for the dogs...
    s = 40
    potentialSSizes = glob.glob(data['alignOutDir'] + '/' + 'miropeats.*ps')
    if len(potentialSSizes) != 0:
        newS = potentialSSizes[0].split('/')[-1]
        newS = newS.split('.')[1]
        s = int(newS)
    data['miropeatSValue'] = s
    
    data['miroOutPS'] = data['alignOutDir'] + '/' + 'miropeats.' + str(s) + '.ps'
    data['miroOutInfo'] = data['alignOutDir'] + '/' + 'miropeats.' + str(s) + '.out'
    if os.path.isfile(data['miroOutPS']) is False:
        data['miroOutPS'] = 'FAILURE'
    
    # to get the gap coordinate file
    get_contig_gaps(data,run=True)
    
    if addBreaks is True:
        data['leftChromFragEnd'] = int(line[7])
        data['leftContigEnd'] = int(line[8])
        data['rightChromFragStart'] = int(line[9])
        data['rightContigStart'] = int(line[10])
###############################################################################
def merge_mr_hits(hitLines):
    newHits = [] 
    for i in hitLines:
        if len(newHits) == 0:
            newHits.append(i)
            continue
        overlapFosmid = genutils.overlap_bed('fos',i[0]-1,i[1],'fos',newHits[-1][0]-1,newHits[-1][1])
        overlapChrom = genutils.overlap_bed('chrom',i[2]-1,i[3],'chrom',newHits[-1][2]-1,newHits[-1][3])
        if overlapFosmid is True and overlapChrom is True:
            if i[1] > newHits[-1][1]:
                 newHits[-1][1] = i[1]
                 newHits[-1][3] = i[3]
        else:
            newHits.append(i)
    return newHits
###############################################################################
def run_3way_align(data):
    figure_align_fragment_to_extract(data)
    get_align_parts(data)
    run_align(data)
    combine_align(data)
    process_3way(data)
    print_fasta_3way(data)    
    get_bp_3way(data)
    print_pretty_alignment(data)
###############################################################################
def repeat_3way_align_after_fragment_update(data,updatedBp = False): # after changed the align parts begin/end manually
    get_align_parts(data)
    run_align(data)
    combine_align(data)
    process_3way(data)
    print_fasta_3way(data)    
    if updatedBp is True:
        get_bp_3way(data)

    print_pretty_alignment(data)
###############################################################################
def populate_previous_align_results(data,line):
    # previous fragments for alignment
    data['genomeFragAlignBegin'] = int(line[11])
    data['genomeFragAlignEnd'] = int(line[12])
    data['contigLeftFragAlignBegin'] = int(line[13])
    data['contigLeftFragAlignEnd'] = int(line[14])
    data['contigRightFragAlignBegin'] = int(line[15])
    data['contigRightFragAlignEnd'] = int(line[16])
    
    # read in the align parts
    get_align_parts(data)
    
    # the align files
    data['leftAlignFileName'] = data['leftSeqFileName'] + '.align'
    data['rightAlignFileName'] = data['rightSeqFileName'] + '.align'
    
    # read in the align file data
    combine_align(data)
    process_3way(data)
    print_fasta_3way(data)
    data['leftBpGenomeFragCoords'] = int(line[17])
    data['leftBpChromCoords'] = int(line[18])
    data['leftBpContigCoord'] = int(line[19])
    
    data['rightBpGenomeFragCoords'] = int(line[20])
    data['rightBpChromCoords'] = int(line[21])
    data['rightBpContigCoord'] = int(line[22])
    
    # why are we recalc instead of reading in?  no reas...
    gTSDs = data['rightBpChromCoords']
    gTSDe = data['leftBpChromCoords']
    gTSDl = gTSDe - gTSDs + 1
    data['gTSDl'] = gTSDl

    insStart = data['leftBpContigCoord'] + 1
    insEnd = data['rightBpContigCoord'] - 1

    s = insEnd-insStart+1
    data['insLen'] = s


    
    
###############################################################################
def populate_all_breakpoint_results(data,line,regDelta,workingDirBase,refGenomeFasta):
    data['refGenomeFasta'] = refGenomeFasta
    populate_data_from_mrfile(data,line,regDelta,workingDirBase,addBreaks=True)
    run_rm(data,run=False) #just to get file names, already ran
    populate_previous_align_results(data,line)
###############################################################################
def figure_align_fragment_to_extract(data):
    # first step, is to setup the contig parts
    extendIn = 60
    extendOut = 60
    
    additDeltaForEnds = 34
    
    # get chrom frag align coordinates

    print data['leftChromFragEnd'],data['rightChromFragStart']
    data['genomeFragAlignBegin'] = data['leftChromFragEnd'] + data['chromFragStart'] - 1 - extendIn + 1
    data['genomeFragAlignEnd'] = data['rightChromFragStart'] + data['chromFragStart'] - 1 + extendIn - 1

    # easy check to make sure there is something -- this is hacky and will lead to manual revisions
    if data['genomeFragAlignEnd'] <= data['genomeFragAlignBegin']:
        t = data['genomeFragAlignBegin']
        data['genomeFragAlignBegin'] = data['genomeFragAlignEnd']
        data['genomeFragAlignEnd'] = t

    # for the contig
    data['contigLeftFragAlignBegin'] = data['leftContigEnd'] - extendOut + 1
    data['contigLeftFragAlignEnd'] = data['leftContigEnd'] + extendOut - 1    
    if data['contigLeftFragAlignBegin'] < 1:
        data['contigLeftFragAlignBegin'] = 1

    if data['contigLeftFragAlignBegin'] == 1 and ((data['contigLeftFragAlignEnd'] - data['contigLeftFragAlignBegin'] +1) <= extendIn):
        data['contigLeftFragAlignEnd'] += additDeltaForEnds


    if data['contigLeftFragAlignEnd'] > data['contigLen']:
        data['contigLeftFragAlignEnd'] = data['contigLen']

    data['contigRightFragAlignBegin'] = data['rightContigStart'] - extendOut + 1
    data['contigRightFragAlignEnd'] = data['rightContigStart'] + extendOut - 1
    if data['contigRightFragAlignEnd'] > data['contigLen']:
        data['contigRightFragAlignEnd'] = data['contigLen']

    if data['contigRightFragAlignEnd'] == data['contigLen'] and ((data['contigRightFragAlignEnd'] - data['contigRightFragAlignBegin'] +1) <= extendIn):
        data['contigRightFragAlignBegin'] -= additDeltaForEnds

    if data['contigRightFragAlignBegin'] < 1:
        data['contigRightFragAlignBegin'] = 1
###############################################################################
def get_align_parts(data):
    data['leftSeqFileName'] =  data['alignOutDir'] + '/' + 'Contig.genomedir.fa.left'
    data['rightSeqFileName'] = data['alignOutDir'] + '/' + 'Contig.genomedir.fa.right'    
    data['genomeFragAlignFileName'] = data['alignOutDir'] + '/' + 'genomeFrag.align.fa'    
    
    contigSeq = genutils.read_fasta_to_string(data['contigSeqFileName'])
    
#    print data['contigLeftFragAlignBegin'],data['contigLeftFragAlignEnd']
#    print data['contigRightFragAlignBegin'],data['contigRightFragAlignEnd']
#    print data['genomeFragAlignBegin'],data['genomeFragAlignEnd']
        
    leftSeq = contigSeq[data['contigLeftFragAlignBegin']-1:data['contigLeftFragAlignEnd']]
    rightSeq = contigSeq[data['contigRightFragAlignBegin']-1:data['contigRightFragAlignEnd']]
    leftSeqStr = genutils.add_breaks_to_line(leftSeq)    
    outFile= open(data['leftSeqFileName'],'w')
    outFile.write('>left\n%s\n' % leftSeqStr)
    outFile.close()
    rightSeqStr = genutils.add_breaks_to_line(rightSeq)    
    outFile= open(data['rightSeqFileName'],'w')
    outFile.write('>right\n%s\n' % rightSeqStr)
    outFile.close()

    # the genome part
    region = data['chromName'] + ':' + str(data['genomeFragAlignBegin']) + '-' + str(data['genomeFragAlignEnd'])
    
    cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + data['genomeFragAlignFileName']
#    print cmd
    genutils.runCMD(cmd)
###############################################################################
def run_align(data):
    data['leftAlignFileName'] = data['leftSeqFileName'] + '.align'
    data['rightAlignFileName'] = data['rightSeqFileName'] + '.align'

    
    cmd = 'stretcher ' + data['leftSeqFileName'] + ' ' + data['genomeFragAlignFileName'] + ' ' + data['leftAlignFileName'] 
    print cmd
    genutils.runCMD(cmd)

    cmd = 'stretcher ' + data['rightSeqFileName'] + ' ' + data['genomeFragAlignFileName'] + ' ' + data['rightAlignFileName'] 
    print cmd
    genutils.runCMD(cmd)
###############################################################################
def read_6_line_block(inFile):
     block = []
     line = inFile.readline()
     if line[0] == '#':
         block.append('END')
         return block
     line = line.rstrip()
     block.append(line) #1
     for i in range(5):
         line = inFile.readline()
         line = line.rstrip()
         block.append(line)
     return block
###############################################################################
def read_align_file(fileName):
    result = {}
    inFile = open(fileName,'r')
    for i in range(15): # take to seq 1 name
        line = inFile.readline()
    # current line is seq 1 name
    line = line.rstrip()
    line = line.split()
    seq1Name = line[-1]
    line = inFile.readline()
    line = line.rstrip()
    line = line.split()
    seq2Name = line[-1]
    result['seq1Name'] =  seq1Name  
    result['seq2Name'] =  seq2Name  
    for i in range(14):
        line = inFile.readline()
    # ok, now we are at the fist block...    
    result['seq1Align'] = ''
    result['seq2Align'] = ''
    
    while True:
        block = read_6_line_block(inFile)
        if block[0] == 'END':
            break

        s1 = block[0].split()[-1]
        s2 = block[2].split()[-1]
    
        if result['seq1Align'] == '':
            result['seq1Align'] = s1
        else:
            result['seq1Align'] += s1
        if result['seq2Align'] == '':
            result['seq2Align'] = s2
        else:
            result['seq2Align'] += s2
    inFile.close()
    return result
###############################################################################
def combine_align(data):
    leftResultsDict = read_align_file(data['leftAlignFileName'])
    rightResultsDict = read_align_file(data['rightAlignFileName'])
    
    genomeSeq = []
    leftSeq = []
    rightSeq = []
    
    leftPtr = 0
    rightPtr = 0
    
    # check of status
    if False:
        print len(leftResultsDict['seq1Align'])
        print len(leftResultsDict['seq2Align'])

        print len(rightResultsDict['seq1Align'])
        print len(rightResultsDict['seq2Align'])
    
    while True:
#        print 'left',leftPtr,'right',rightPtr
        # see if done with both.
        # for now, we are assuming left and right start align start inside of genome, so we are accounting for all
        if leftPtr >= len(leftResultsDict['seq1Align']) and rightPtr >= len(rightResultsDict['seq1Align']):
            break
        
        if leftPtr < len(leftResultsDict['seq1Align']):
            leftS1 = leftResultsDict['seq1Align'][leftPtr]
            leftS2 = leftResultsDict['seq2Align'][leftPtr]
        else:
            leftS1 = '-'
            leftS2 = '-'
            
        
        if rightPtr < len(rightResultsDict['seq1Align']):
            rightS1 = rightResultsDict['seq1Align'][rightPtr]
            rightS2 = rightResultsDict['seq2Align'][rightPtr]
        else:
            rightS1 = '-'
            rightS2 = '-'
        
        # same genome seq
        if leftS2 == rightS2:
            genomeSeq.append(leftS2)
            leftSeq.append(leftS1)
            rightSeq.append(rightS1)
            leftPtr += 1
            rightPtr += 1
            continue
        
        # gap in left genome, not in left contig
        if leftS2 == '-' and leftS1 != '-':
            genomeSeq.append(leftS2)
            leftSeq.append(leftS1)
            rightSeq.append('-')
            leftPtr += 1
            continue
        # gap in right genome, not in right contig    
        if rightS2 == '-' and rightS1 != '-':
            genomeSeq.append(rightS2)
            leftSeq.append('-')
            rightSeq.append(rightS1)
            rightPtr += 1
            continue
        print 'how did I get here?'
        sys.exit()
#    print len(genomeSeq),len(leftSeq),len(rightSeq)
    
    data['left3way'] = leftSeq
    data['right3way'] = rightSeq
    data['genome3way'] = genomeSeq
###############################################################################
# goal is to setup and process the 3 way aignment data
def process_3way(data):
    data['left3wayPos'] = []
    prev = 0
    for i in range(len(data['left3way'])):
        if data['left3way'][i] == '-':
            data['left3wayPos'].append(prev)
        else:
            prev += 1
            data['left3wayPos'].append(prev)
            
    data['right3wayPos'] = []
    prev = 0
    for i in range(len(data['right3way'])):
        if data['right3way'][i] == '-':
            data['right3wayPos'].append(prev)
        else:
            prev += 1
            data['right3wayPos'].append(prev)

    data['genome3wayPos'] = []
    prev = 0
    for i in range(len(data['genome3way'])):
        if data['genome3way'][i] == '-':
            data['genome3wayPos'].append(prev)
        else:
            prev += 1
            data['genome3wayPos'].append(prev)
    # no, do assignment
    # 1 == left matches
    # 2 == right matches
    # * all 3 matches
    # ' ' both left and right are gaps
    # 'N' not gap, but not match either        
    
    data['3wayParse'] = []
    for i in range(len(data['genome3way'])):
        l = data['left3way'][i]
        r = data['right3way'][i]
        g = data['genome3way'][i]
        if g == l and g == r:
            data['3wayParse'].append('*')
        elif g == l and g != r :  
            data['3wayParse'].append('1')
        elif g != l and g == r :
            data['3wayParse'].append('2')
        elif l == '-' and r == '-':
            data['3wayParse'].append(' ')
        else:
            data['3wayParse'].append('N')

###############################################################################
def print_fasta_3way(data):
    # print out as a fasta for later
    data['3wayAlignFileFAName'] = data['alignOutDir'] + '/' + data['siteID']  + '.3wayalign.fasta'    


    outFile = open(data['3wayAlignFileFAName'],'w')
    nl = ''.join(data['left3way'])    
    outFile.write('>left\n%s\n' % nl)
    nl = ''.join(data['genome3way'])    
    outFile.write('>chrom\n%s\n' % nl)
    nl = ''.join(data['right3way'])    
    outFile.write('>right\n%s\n' % nl)
    # print out as a fasta
    outFile.close()
###############################################################################
# go through the 3 way data and get the breakpoints
def get_bp_3way_lasthit(data):
    # get last 1 before we hit 2'sss
    # get first_p
    first1 = 0
    for first1 in range(len(data['3wayParse'])):
        if data['3wayParse'][first1] == '1':
            break
    # now, continue on till we get the first 2
    for first2 in range(first1,len(data['3wayParse'])):
        if data['3wayParse'][first2] == '2':
            break
    # now, go backwards from there to find the last1
    for last1 in range(first2,0,-1):
        if data['3wayParse'][last1]  in['1','*']:
            break
    print 'have last1',last1
    leftBp = last1

   # get first 2, 
    last2 = len(data['3wayParse']) -1
    for last2 in range(len(data['3wayParse'])-1,0,-1):
        if data['3wayParse'][last2] == '2':
            break
    # continue on until we get a 1
    for last1 in range(last2,0,-1):
        if data['3wayParse'][last1] == '1':
            break
    # go forward till we have a 2
    for first2 in range(last1,len(data['3wayParse'])):
        if data['3wayParse'][first2]  in ['*','2']:
            break
    
    print 'have first2,',first2
    rightBp = first2
        
    
    data['leftBpCol'] = leftBp
    data['rightBpCol'] = rightBp
    
    data['leftBpContigCoord'] = data['left3wayPos'][data['leftBpCol']] + data['contigLeftFragAlignBegin'] - 1
    data['rightBpContigCoord'] = data['right3wayPos'][data['rightBpCol']] + data['contigRightFragAlignBegin'] - 1
    
    data['leftBpGenomeFragCoords'] = data['genome3wayPos'][data['leftBpCol']]
    data['rightBpGenomeFragCoords'] = data['genome3wayPos'][data['rightBpCol']]


    print 'leftbp',data['leftBpContigCoord'],data['leftBpGenomeFragCoords']
    print 'rightbp',data['rightBpContigCoord'],data['rightBpGenomeFragCoords']

    data['leftBpChromCoords'] = data['leftBpGenomeFragCoords'] + data['genomeFragAlignBegin'] - 1
    data['rightBpChromCoords'] = data['rightBpGenomeFragCoords'] + data['genomeFragAlignBegin'] - 1
    
#    for i in range(leftBp-10,rightBp+10):
#         print i,data['left3way'][i],data['genome3way'][i],data['right3way'][i],data['3wayParse'][i]
###############################################################################
###############################################################################
# go through the 3 way data and get the breakpoints
# This version is based on idea of maximum score
# +1 for match, -1 for mismatch, then take the max
#updates on 20 August for dealing with N "GAP" sequence in alignments
def get_bp_3way(data):
    # get last 1 before we hit 2'sss
    # get first_p
    matchScore = 1
    mismatchScore = -3
    bothScore = 0
#    neitherScore = 0
    neitherScore = -1

    leftScores = []    
    foundFirst = False
    first1 = 0
    for first1 in range(len(data['3wayParse'])):
        k = data['3wayParse'][first1]
        if k == '1' or k == '*':
            foundFirst = True
            if len(leftScores) == 0:
                leftScores.append(matchScore)
            else:
                leftScores.append(leftScores[-1] + matchScore)
        else:
            if k == 'N':
                s = neitherScore
            else:
                s = mismatchScore
            if foundFirst is False:

                leftScores.append(s)
            else:
                leftScores.append(max(0,leftScores[-1] + s))


    rightScores = []    
    foundLast = False
    last1 = 0

    for last1 in range(len(data['3wayParse'])-1,-1,-1):
        k = data['3wayParse'][last1]
        if k == '2' or k == '*':
            foundLast = True
            if len(rightScores) == 0:
                rightScores.append(matchScore)
            else:
                rightScores.append(rightScores[-1] + matchScore)
        else:
            if k == 'N':
                s = neitherScore
            else:
                s = mismatchScore
            
            if foundLast is False:
                rightScores.append(s)
            else:
                rightScores.append(max(0,rightScores[-1] + s))
    rightScores.reverse() #since we filled it in backwards
        
    data['leftBpCol'] =  genutils.get_max_index(leftScores)
    data['rightBpCol'] = genutils.get_max_index(rightScores)
    
 #   print 'max for left is',data['leftBpCol']
 #   print 'max for right is',data['rightBpCol']
 #   for i in range(len(data['3wayParse'])):
 #       print i,data['3wayParse'][i],leftScores[i],rightScores[i]
 
 
        
#    sys.exit()
    


    data['leftBpContigCoord'] = data['left3wayPos'][data['leftBpCol']] + data['contigLeftFragAlignBegin'] - 1
    data['rightBpContigCoord'] = data['right3wayPos'][data['rightBpCol']] + data['contigRightFragAlignBegin'] - 1
    
    data['leftBpGenomeFragCoords'] = data['genome3wayPos'][data['leftBpCol']]
    data['rightBpGenomeFragCoords'] = data['genome3wayPos'][data['rightBpCol']]


    print 'leftbp',data['leftBpContigCoord'],data['leftBpGenomeFragCoords']
    print 'rightbp',data['rightBpContigCoord'],data['rightBpGenomeFragCoords']

    data['leftBpChromCoords'] = data['leftBpGenomeFragCoords'] + data['genomeFragAlignBegin'] - 1
    data['rightBpChromCoords'] = data['rightBpGenomeFragCoords'] + data['genomeFragAlignBegin'] - 1
    
   
    

#    for i in range(leftBp-10,rightBp+10):
#         print i,data['left3way'][i],data['genome3way'][i],data['right3way'][i],data['3wayParse'][i]
###############################################################################
# print out alignment as pretty sequence
def print_pretty_alignment(data):    
    #left end is blue, right start is red
    data['3wayAlignFilePrettyName'] = data['genomeFragFileName'] + '.3wayalign.pretty'    
    data['3wayAlignFilePrettyNamePS'] = data['alignOutDir'] + '/' + data['siteID'] + '.3wayalign.pretty.ps'
    data['3wayAlignFilePrettyNamePDF'] = data['alignOutDir'] + '/' + data['siteID'] + '.3wayalign.pretty.pdf'


    outFile = open(data['3wayAlignFilePrettyName'],'w')
    outFile.write('Site ID: %s\n' % (data['siteID']))
    outFile.write('%s\t%s\n' % (data['contigName'],data['contigDir']))
    outFile.write('%s:%i-%i\n' % (data['chromName'],data['genomeFragAlignBegin'],data['genomeFragAlignEnd']))
    # left BP in chromFrag and Contig
    outFile.write('~color{0 0 1}end left match~color{default} chromFrag %i Contig %i\n' % (data['leftBpGenomeFragCoords'],data['leftBpContigCoord']))
    # right BP in chromFrag and Contig
    outFile.write('~color{1 0 0}start right match~color{default} chromFrag %i Contig %i\n' % (data['rightBpGenomeFragCoords'],data['rightBpContigCoord']))


    #go through and add in the colors
    # do the colors individually
    print 'ready to start'
    print data['leftBpContigCoord'],data['leftBpGenomeFragCoords']
    print data['rightBpContigCoord'],data['rightBpGenomeFragCoords']
    
    

    
    
    for i in range(0,len(data['genome3way'])):
        if data['left3wayPos'][i] == (data['leftBpContigCoord'] - data['contigLeftFragAlignBegin'] + 1):
            if data['left3way'][i] == '-':
                print 'left is -'
            else:
                data['left3way'][i] = '~color{0 0 1}' + data['left3way'][i] + '~color{default}'
                data['3wayParse'][i] = '~color{0 0 1}' + data['3wayParse'][i] + '~color{default}'
                print 'LEFT contig',i,data['left3wayPos'][i]

        if data['genome3wayPos'][i] == data['leftBpGenomeFragCoords'] and data['genome3way'][i] != '-':
            data['genome3way'][i] = '~color{0 0 1}' + data['genome3way'][i] + '~color{default}'
            print 'LEFT GENOME',i,data['genome3wayPos'][i]

        
        if data['right3wayPos'][i] == (data['rightBpContigCoord'] - data['contigRightFragAlignBegin'] +1):
            if data['right3way'][i] == '-':
                print i,'right is -'
            else:
                data['right3way'][i] = '~color{1 0 0}' + data['right3way'][i] + '~color{default}'
                data['3wayParse'][i] = '~color{1 0 0}' + data['3wayParse'][i] + '~color{default}'
                print i,'right contig',i,data['right3wayPos'][i]
        if data['genome3wayPos'][i] == data['rightBpGenomeFragCoords'] and data['genome3way'][i] != '-' :
            data['genome3way'][i] = '~color{1 0 0}' + data['genome3way'][i] + '~color{default}'
            print i,'RIGHT GENOME',i,data['genome3wayPos'][i]
        
    leftName =  'left   '
    rightName = 'right  '
    chromName = 'chrom  '
    passeName = '       '
    outFile.write('\n\n')
    # do it in runs of 50
    width = 70
    sliceS = 0
    sliceE = sliceS + width
    while True:
        if sliceS >= len(data['genome3way']):
            break
        if sliceE > len(data['genome3way']):
            sliceE = len(data['genome3way'])
        l = data['left3way'][sliceS:sliceE]
        g = data['genome3way'][sliceS:sliceE]
        r = data['right3way'][sliceS:sliceE]
        p = data['3wayParse'][sliceS:sliceE]
        
        l = leftName + ''.join(l)
        g = chromName + ''.join(g)
        r = rightName + ''.join(r)
        p = passeName + ''.join(p)
        
        outFile.write('%s\n%s\n%s\n%s\n\n' % (l,g,r,p))
        sliceS = sliceE
        sliceE = sliceS + width
    outFile.close()    

    print 'Clean up PS and PDF'
    cmd = 'rm ' + data['3wayAlignFilePrettyNamePS']
    print cmd
    genutils.runCMDNoFail(cmd)
    cmd = 'rm ' + data['3wayAlignFilePrettyNamePDF']
    print cmd
    genutils.runCMDNoFail(cmd)
    

    cmd = 'enscript %s -o %s -e~ -B -2r' % (data['3wayAlignFilePrettyName'],data['3wayAlignFilePrettyNamePS'])
    print cmd
    genutils.runCMD(cmd)
    
    cmd = 'ps2pdf ' + data['3wayAlignFilePrettyNamePS'] + ' ' + data['3wayAlignFilePrettyNamePDF']
    print cmd
    genutils.runCMD(cmd)
###############################################################################
def write_insertion_genome_dir(data,outRef):
    insSeq = data['contigSeqGenomeDir'][data['leftBpContigCoord']:data['rightBpContigCoord']-1] 
#    print len(insSeq)
#    print insSeq
#    print data['siteID']
    outS = genutils.add_breaks_to_line(insSeq)
    outRef.write('>%s\n%s\n' % (data['siteID'],outS))    
###############################################################################
def write_insertion_repeat_dir(data,outRepeat):
    insSeq = data['contigSeqGenomeDir'][data['leftBpContigCoord']:data['rightBpContigCoord']-1] 
    rmLines = read_rm_file(data['contigSeqFileNameRM'])
    # figure if sine is on + or -
    rmDir = figure_sine_strand(rmLines)
    if rmDir == '-':
        insSeq = genutils.revcomp(insSeq)
    
    outS = genutils.add_breaks_to_line(insSeq)
    outRepeat.write('>%s\n%s\n' % (data['siteID'],outS))    
###############################################################################
def read_rm_file(rmFileName):
    rmLines = []
    inFile = open(rmFileName,'r')
    for line in inFile:
        if line == '\n':
            continue
        line = line.rstrip()
        line = line.split()
        if line[0] == 'There':
            return []
        if line[0] == 'SW':
            continue
        if line[0] == 'score':
            continue
        rmLines.append(line)
    inFile.close()
    return rmLines
###############################################################################
def figure_sine_strand(rmLines):
    numPlus = 0
    numMinus = 0
    for i in rmLines:
        if i[10] == 'SINE/Alu':
            if i[8] == '+':
                numPlus += 1
            else:
                numMinus += 1
    if numPlus > numMinus:
        return '+'
    else:
        return '-'
###############################################################################
def make_alternative_seqs(data,bpOutTable,allelesBaseDir,fragmentExtension):
    alleleDir = allelesBaseDir + data['siteID']
    if os.path.isdir(alleleDir) is False:
        cmd = 'mkdir ' + alleleDir
        print cmd
        genutils.runCMD(cmd)
    
    alleleDir += '/'
    genomeLeftFa = alleleDir + 'genomeLeft.fa'
    genomeRightFa = alleleDir + 'genomeRight.fa'
    genomeWholeFa = alleleDir + 'genomeWhole.fa'    
    alleleFa = alleleDir + 'alleles.fa'
    data['alleleFa'] = alleleFa
    data['alleleDir'] = alleleDir
    gTSDs = data['rightBpChromCoords']
    gTSDe = data['leftBpChromCoords']
    gTSDl = gTSDe - gTSDs + 1
    data['gTSDl'] = gTSDl

    if gTSDl <= -1:  # deletion in chromosome
        print 'deletion of %i in genome' % gTSDl
        leftChromBp = data['leftBpChromCoords']
        leftChromStart = leftChromBp - fragmentExtension + 1
        rightChromBp = data['rightBpChromCoords']
        rightChromEnd = rightChromBp + fragmentExtension - 1       
        data['insSite'] = leftChromBp
                
#        print leftChromBp,leftChromStart,rightChromBp,rightChromEnd    
        region = data['chromName'] + ':' + str(leftChromStart) + '-' + str(leftChromBp)
        cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + genomeLeftFa
        genutils.runCMD(cmd)
        genomeLeftSeq = genutils.read_fasta_to_string(genomeLeftFa)
        genomeLeftSeq = genomeLeftSeq.upper()
        region = data['chromName'] + ':' + str(rightChromBp) + '-' + str(rightChromEnd)
        cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + genomeRightFa
        genutils.runCMD(cmd)
        genomeRightSeq = genutils.read_fasta_to_string(genomeRightFa)
        genomeRightSeq = genomeRightSeq.upper()
        # get the chrom sequence
        region = data['chromName'] + ':' + str(leftChromStart) + '-' + str(rightChromEnd)
        cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + genomeWholeFa
        genutils.runCMD(cmd)
        genomeWholeSeq = genutils.read_fasta_to_string(genomeWholeFa)
        genomeWholeSeq = genomeWholeSeq.upper()
        
        data['mapFragStart'] = leftChromStart
        data['mapFragEnd'] = rightChromEnd
        # since that BP is in contig
        contigStart = data['leftBpContigCoord'] # already last bp
        contigEnd = data['rightBpContigCoord']         
        contigSeq = data['contigSeqGenomeDir'][contigStart:contigEnd-1]

        # print out genome        
        outFile = open(alleleFa,'w')
        outFile.write('>%s\n' % (data['siteID']+'_genome'))
        gSeq = genomeWholeSeq
        gSeq = genutils.add_breaks_to_line(gSeq)
        outFile.write('%s\n' % gSeq)
        outFile.write('>%s\n' % (data['siteID']+'_insertion'))
        iSeq = genomeLeftSeq + contigSeq + genomeRightSeq
        iSeq = genutils.add_breaks_to_line(iSeq)
        outFile.write('%s\n' % iSeq)
        outFile.close()
#        print 'left',genomeLeftSeq
#        print 'right',genomeRightSeq
#        print 'contig',contigSeq
#        print len(contigSeq)
    elif  gTSDl >= 1:  # has TSD
#        print 'has TSD len %i' % gTSDl
        # note that they cross
        leftChromBp = data['rightBpChromCoords']
        leftChromStart = leftChromBp - fragmentExtension + 1
        rightChromBp = data['leftBpChromCoords']
        rightChromEnd = rightChromBp + fragmentExtension - 1
        
        data['insSite'] = leftChromBp
#        print leftChromBp,leftChromStart,rightChromBp,rightChromEnd    
        region = data['chromName'] + ':' + str(leftChromStart) + '-' + str(leftChromBp)
        cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + genomeLeftFa
        genutils.runCMD(cmd)
        genomeLeftSeq = genutils.read_fasta_to_string(genomeLeftFa)
        genomeLeftSeq = genomeLeftSeq.upper()
        region = data['chromName'] + ':' + str(rightChromBp) + '-' + str(rightChromEnd)
        cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + genomeRightFa
        genutils.runCMD(cmd)
        genomeRightSeq = genutils.read_fasta_to_string(genomeRightFa)
        genomeRightSeq = genomeRightSeq.upper()
        # get the chrom sequence
        region = data['chromName'] + ':' + str(leftChromStart) + '-' + str(rightChromEnd)
        cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + genomeWholeFa
        genutils.runCMD(cmd)
        
        data['mapFragStart'] = leftChromStart
        data['mapFragEnd'] = rightChromEnd
        genomeWholeSeq = genutils.read_fasta_to_string(genomeWholeFa)
        genomeWholeSeq = genomeWholeSeq.upper()
        
        # since that BP is in contig
        contigStart = data['leftBpContigCoord']  - gTSDl + 1  # to get over to right size, include the TSD
        contigEnd = data['rightBpContigCoord']  + gTSDl  - 1 # to get over to the right size, include the TSD        
        contigSeq = data['contigSeqGenomeDir'][contigStart:contigEnd-1]

        # print out genome        
        outFile = open(alleleFa,'w')
        outFile.write('>%s\n' % (data['siteID']+'_genome'))
        gSeq = genomeWholeSeq
        gSeq = genutils.add_breaks_to_line(gSeq)
        outFile.write('%s\n' % gSeq)
        outFile.write('>%s\n' % (data['siteID']+'_insertion'))
        iSeq = genomeLeftSeq + contigSeq + genomeRightSeq
        iSeq = genutils.add_breaks_to_line(iSeq)
        outFile.write('%s\n' % iSeq)
        outFile.close()
#        print 'left',genomeLeftSeq
#        print 'right',genomeRightSeq
#        print 'contig',contigSeq
#        print len(contigSeq)
    elif  gTSDl == 0:  # has no TSD 
#        print 'has TSD len %i' % gTSDl
        # note that they cross
        leftChromBp = data['leftBpChromCoords']
        leftChromStart = leftChromBp - fragmentExtension + 1
        rightChromBp = data['rightBpChromCoords']
        rightChromEnd = rightChromBp + fragmentExtension - 1        
        data['insSite'] = leftChromBp       
#        print leftChromBp,leftChromStart,rightChromBp,rightChromEnd    
        region = data['chromName'] + ':' + str(leftChromStart) + '-' + str(leftChromBp)
        cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + genomeLeftFa
        genutils.runCMD(cmd)
        genomeLeftSeq = genutils.read_fasta_to_string(genomeLeftFa)
        genomeLeftSeq = genomeLeftSeq.upper()
        region = data['chromName'] + ':' + str(rightChromBp) + '-' + str(rightChromEnd)
        cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + genomeRightFa
        genutils.runCMD(cmd)
        genomeRightSeq = genutils.read_fasta_to_string(genomeRightFa)
        genomeRightSeq = genomeRightSeq.upper()
        # get the chrom sequence
        region = data['chromName'] + ':' + str(leftChromStart) + '-' + str(rightChromEnd)
        cmd = 'samtools faidx ' + data['refGenomeFasta'] + ' ' + region + ' > ' + genomeWholeFa
        genutils.runCMD(cmd)        
        data['mapFragStart'] = leftChromStart
        data['mapFragEnd'] = rightChromEnd        
        genomeWholeSeq = genutils.read_fasta_to_string(genomeWholeFa)
        genomeWholeSeq = genomeWholeSeq.upper()
        
        # since that BP is in contig
        contigStart = data['leftBpContigCoord']  
        contigEnd = data['rightBpContigCoord']          
        contigSeq = data['contigSeqGenomeDir'][contigStart:contigEnd-1]

        # print out genome        
        outFile = open(alleleFa,'w')
        outFile.write('>%s\n' % (data['siteID']+'_genome'))
        gSeq = genomeWholeSeq
        gSeq = genutils.add_breaks_to_line(gSeq)
        outFile.write('%s\n' % gSeq)
        outFile.write('>%s\n' % (data['siteID']+'_insertion'))
        iSeq = genomeLeftSeq + contigSeq + genomeRightSeq
        iSeq = genutils.add_breaks_to_line(iSeq)
        outFile.write('%s\n' % iSeq)
        outFile.close()
#        print 'left',genomeLeftSeq
#        print 'right',genomeRightSeq
#        print 'contig',contigSeq
#        print len(contigSeq)
    else:
        print 'What TSD size?'
        print gTSDl
        sys.exit()
    # make out file
    nl = [data['siteID'],data['chromName'],data['insSite'],gTSDl,data['mapFragStart'],data['mapFragEnd'] ]
    nl = [str(i) for i in nl]
    nl = '\t'.join(nl) + '\n'
    bpOutTable.write(nl)
    bwa_index_alleles(data)
###############################################################################
def bwa_index_alleles(data):
    cmd = 'bwa-0.5.9 index %s' % (data['alleleFa'])
#    print cmd
    genutils.runCMD(cmd)
###############################################################################
def process_gap_file(fn):
    gapLines = []
    inFile = open(fn,'r')
    for line in inFile:
        line = line.rstrip()
        line = line.split()
        s = int(line[1])
        e = int(line[2])
        gapLines.append([s,e])
    inFile.close()
    return gapLines
###############################################################################
def regulate_align_frags(data):
    if data['contigLeftFragAlignBegin'] < 1:
        data['contigLeftFragAlignBegin' ] = 1
    if data['contigLeftFragAlignEnd'] > data['contigLen'] :
        data['contigLeftFragAlignEnd'] = data['contigLen']

    if data['contigRightFragAlignBegin'] < 1:
        data['contigRightFragAlignBegin' ] = 1
    if data['contigRightFragAlignEnd'] > data['contigLen'] :
        data['contigRightFragAlignEnd'] = data['contigLen']
###############################################################################
