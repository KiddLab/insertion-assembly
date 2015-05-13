import sys
import pickle
from optparse import  OptionParser
scriptDir = '/home/jmkidd/kidd-lab/jmkidd-projects/people-projects/jwilds-projects/RetroSeq-HGDP/retroseq-alu-additional/results/assembly/hg18/scripts/'


###############################################################################
#changes based on PS_template_April2010 value
def display_edit_setup(line):
     if '%%DocumentFonts: Helvetica-Bold Helvetica' in line:
         line = '%%DocumentFonts: Arial-Bold Arial\n'
     
     if '%%BoundingBox:' in line:
         line = '%%BoundingBox: 0 0 612 792\n'
     if '/repwidth' in line:
         line = '/repwidth 0.25 cm def\n'
     if '/linkwidth' in line:
         line = '/linkwidth 0 def\n'
     if '/gapwidth' in line:
         line = '/gapwidth 0 def\n'
     if '/leftmargin' in line:
         line = '/leftmargin 2 cm def\n'
     if '/bottommargin' in line:
         line = '/bottommargin 2 cm def\n'
     if '/pageheight' in line:
         line = '/pageheight 26 cm def\n'

     if '/Helvetica-Bold findfont' in line:
         line = '/Arial-Bold findfont 12 scalefont setfont\n\n'

          
     if '/Helvetica findfont' in line:
         line = '/Arial findfont 9 scalefont setfont\n'
     if '0 pageheight 2 cm sub moveto' in line:
         line = '0 pageheight 1.5 cm sub moveto\n'
     if '/graphicmargin 17.5' in line:
         line = '/graphicmargin 17.5 cm def\n'
         
     if '(Longest Sequence ' in line:
         line = '\n'    
     if '( Threshold Score ' in line:
         line = '\n'    

         
     return line
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
def repeat_class_to_name(r):
    if 'SINE' in r:
        return 'SINE'
    if 'ARTEFACT' in r:
        return 'ARTEFACT'
    if 'DNA' in r:
        return 'DNA'
    if 'LINE' in r:
        return 'LINE'
    if 'Low_complexity' in r:
        return 'Low_complexity'
    if 'LTR' in r:
        return 'LTR'
    if 'Other' in r:
        return 'Other'
    if 'rRNA' in r:
        return 'rRNA'
    if 'scRNA' in r:
        return 'scRNA'
    if 'snRNA' in r:
        return 'snRNA'
    if 'srpRNA' in r:
        return 'srpRNA'
    if 'tRNA' in r:
        return 'tRNA'
    if 'RNA' in r:
        return 'RNA'
    if 'Satellite' in r:
        return 'Satellite'
    if 'Simple_repeat' in r:
        return 'Simple_repeat'
    if 'Unknown' in r:
        return 'Unknown'
    if 'Retroposon' in r:
        return 'SINE'    
    if 'RC/Helitron' in r:    
        return 'Other'
    if 'Helitron' in r:    
        return 'Other'

    print 'repeat class unknown for',r
    sys.exit()
###############################################################################
def sort_and_merge_repeats(repeats):
    repeats.sort()
    newRep = []
    for r in repeats:
        s = r[0]
        e = r[1]
        orientation = r[2]
        repClass = r[3]
        if len(newRep) == 0:
            newRep.append(r)
        else:
            lr = newRep[-1]
            ls = lr[0]
            le = lr[1]
            lorient = lr[2]
            lrepClass = lr[3]
            # overlap, need to extend
            if le > s and lrepClass == repClass and lorient == orientation:
                n = [ls,e,orientation,repClass]
                newRep[-1] = n
            else:
                newRep.append(r)
    return newRep    
###############################################################################
def process_repeat_file(fn):
    repeatLines = read_rm_file(fn)
    repeats = []
    for R in repeatLines:
        s = int(R[5])
        e = int(R[6])
        orientation = R[8]
        repClass = R[10]
        if orientation == 'C':
            orientation = '-'
        repClass = repeat_class_to_name(repClass)
        reps = [s,e,orientation,repClass]
        repeats.append(reps)
    repeats = sort_and_merge_repeats(repeats)
    return repeats
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


USAGE = """
python annotate-miropeats-3align.py    --pickle <pickle file name>  


"""
###############################################################################
parser = OptionParser(USAGE)
parser.add_option('--pickle',dest='dataPickleName',help='input data pickle name')

(options,args)=parser.parse_args()
if options.dataPickleName is None:
    parser.error('data pickle file name name not given')
###############################################################################



# Info files
PSTemplate = scriptDir +'/PS_template_July2008_withlines.txt' 
#### setup the repeat to color dictionary#####
color = {}
color['Other'] = 'Black'
color['Simple Repeat'] = 'DkGray'
color['Low Complexity'] = 'LtGray'
color['DNA'] = 'Pink'
color['LTR'] = 'Orange'
color['LINE'] = 'Green'
color['SINE'] = 'Purple'
#height of ref and fos sequence, in "Y" coordinates
#fos is what is on the bottom
ref_line=0.6
fos_line=0.1

# for direction...
#arrowEndBp = 1000
arrowEndBp = 100


###############################################################################
print 'reading in data from',options.dataPickleName
data = read_pickle_dictionary(options.dataPickleName)



###############################################################################


miropeatsInFile = data['miroOutPS']
topRepeatMaskFile = data['genomeFragFileNameRM']
bottomRepeatMaskFile = data['contigSeqFileNameRM']


miropeatsOutFile = miropeatsInFile + '.breakpoints.annotated.ps'

# these are for coloring breakpoints... we will skip this for now...
chrm_breaks = {}
chrm_breaks['not do'] = 1
clone_breaks = {}
clone_breaks['not do'] = 1


inFile = open(miropeatsInFile,'r')
outFile = open(miropeatsOutFile,'w')


# read in to tagends function, then add in the template
while True:
    line = inFile.readline()
    if line[0:8] == '/tagends':
        break
    else:
        line = display_edit_setup(line)
        outFile.write(line)
#here, line is the begin of the tagends function

# here is a good place to print out other info
print 'Adding PS template info.....'
inTEMP = open(PSTemplate,'r')
for t in inTEMP:
    outFile.write(t)
inTEMP.close()
print 'added the template info!'



###############################################################################
# This is for printing out the pretty info
outFile.write('/Arial-Bold findfont 12 scalefont setfont\n')
outFile.write('0 pageheight 0.6 cm sub moveto\n')
outFile.write('(siteID: %s %s %s) show\n' % (data['siteID'],data['contigName'],data['contigDir']))
outFile.write('0 pageheight 1.1 cm sub moveto\n')
outFile.write('(%s:%i-%i) show\n' % (data['chromName'],data['genomeFragAlignBegin'],data['genomeFragAlignEnd']))
outFile.write('0 pageheight 1.6 cm sub moveto\n')
outFile.write('(TSD size: %s INS size: %s) show\n' % (data['gTSDl'],data['insLen']) )
outFile.write('/Courier findfont 10 scalefont setfont\n')
y = 2.1
outFile.write('0 pageheight %f cm sub moveto\n' % y)

leftName =  'left   '
rightName = 'right  '
chromName = 'chrom  '
passeName = '       '
# do it in runs of 50
width = 70
sliceS = 0
sliceE = sliceS + width
print 'ADDING ALIGNMENT TEXT'
textLineDelta = 0.3
while True:
    if sliceS >= len(data['genome3way']):
        break
    if sliceE > len(data['genome3way']):
        sliceE = len(data['genome3way'])
    l = data['left3way'][sliceS:sliceE]
    g = data['genome3way'][sliceS:sliceE]
    r = data['right3way'][sliceS:sliceE]
    p = data['3wayParse'][sliceS:sliceE]
    
    y += textLineDelta
    outFile.write('0 pageheight %f cm sub moveto\n' % y)
    outFile.write('(%s) show\n' % (leftName))
    for  i in l:
        if i[0] == '~':
            if '~color{0 0 1}' in i:
                outFile.write('Blue\n')
                outFile.write('(%s) show\n' % (i[13]))
                outFile.write('Black\n')
            elif '~color{1 0 0}' in i:
                outFile.write('Red\n')
                outFile.write('(%s) show\n' % (i[13]))
                outFile.write('Black\n')
            else:
                print 'What color?'
                print i
                sys.exit()
        else:
            outFile.write('(%s) show\n' % (i))
    y += textLineDelta

    outFile.write('0 pageheight %f cm sub moveto\n' % y)
    outFile.write('(%s) show\n' % (chromName))
    for  i in g:
        if i[0] == '~':
            if i.count('~') > 2:
                # we have multiple colors for the same pos
                outFile.write('(%s) show\n' % (i[26]))
                continue

            if '~color{0 0 1}' in i:
                outFile.write('Blue\n')
                outFile.write('(%s) show\n' % (i[13]))
                outFile.write('Black\n')
            elif '~color{1 0 0}' in i:
                outFile.write('Red\n')
                outFile.write('(%s) show\n' % (i[13]))
                outFile.write('Black\n')
            else:
                print 'What color?'
                print i
                sys.exit()
        else:
            outFile.write('(%s) show\n' % (i))
    y += textLineDelta

    outFile.write('0 pageheight %f cm sub moveto\n' % y)
    outFile.write('(%s) show\n' % (rightName))
    for  i in r:
        if i[0] == '~':
            if '~color{0 0 1}' in i:
                outFile.write('Blue\n')
                outFile.write('(%s) show\n' % (i[13]))
                outFile.write('Black\n')
            elif '~color{1 0 0}' in i:
                outFile.write('Red\n')
                outFile.write('(%s) show\n' % (i[13]))
                outFile.write('Black\n')
            else:
                print 'What color?'
                print i
                sys.exit()
        else:
            outFile.write('(%s) show\n' % (i))
    y +=textLineDelta

    outFile.write('0 pageheight %f cm sub moveto\n' % y)
    outFile.write('(%s) show\n' % (passeName))
    for  i in p:
        if i[0] == '~':            
            if i.count('~') > 2:
                # we have multiple colors for the same pos
                outFile.write('(%s) show\n' % (i[26]))
                continue
            if '~color{0 0 1}' in i:
                outFile.write('Blue\n')
                outFile.write('(%s) show\n' % (i[13]))
                outFile.write('Black\n')
            elif '~color{1 0 0}' in i:
                outFile.write('Red\n')
                outFile.write('(%s) show\n' % (i[13]))
                outFile.write('Black\n')
            else:
                print 'What color?'
                print i
                sys.exit()
        else:
            outFile.write('(%s) show\n' % (i))
    y += textLineDelta
    sliceS = sliceE
    sliceE = sliceS + width


outFile.write('/Arial findfont 12 scalefont setfont\n')    
###############################################################################






# print out rest of info down to where coords will go
outFile.write(line)
started_pc = 0
more  = 1
while more == 1:
    line = inFile.readline()
    
    if 'printcontig' in line:
        line = line.rstrip()
        started_pc = 1
        f = line.split()
        # determine which one is ref line, it has '-' in it

        
        if '-' in f[1] :
            f[0] = ref_line
            f[1] = f[1].split(':')[0] + ')'
            
        else:
            f[0] = fos_line
        
        if '-' in f[6] :
            f[5] = ref_line
            f[6] = f[6].split(':')[0] + ')'

        else:
            f[5] = fos_line
        p_cmd = 'printcontig'

        # The below is if we are going to color potential chromosome breakpoints that we have called based on the
        # miropeats output directly
        #chrom first
        if 'chr' in f[1]:
            if f[2] in chrm_breaks and f[7] in clone_breaks:
                p_cmd = 'printcontig_left'
            if f[3] in chrm_breaks and f[8] in clone_breaks:
                p_cmd = 'printcontig_right'
            if f[3] in chrm_breaks and f[8] in clone_breaks and f[2] in chrm_breaks and f[7] in clone_breaks:
                p_cmd = 'printcontig_both'
        else:
            if f[7] in chrm_breaks and f[2] in clone_breaks:
                p_cmd = 'printcontig_left'
            if f[8] in chrm_breaks and f[3] in clone_breaks:
                p_cmd = 'printcontig_right'
            if f[3] in chrm_breaks and f[8] in clone_breaks and f[2] in chrm_breaks and f[7] in clone_breaks:
                p_cmd = 'printcontig_both'
        
        print p_cmd
        f[10] = p_cmd
        f = [str(j) for j in f]
        line = ' '.join(f) + '\n'
    if (('printcontig' in line) is False) and (started_pc == 1):  # means that we are done with the printcontig commands
        more = 0
    if more == 1:
        outFile.write(line)
        
#OK now time for repeats info                       
print 'Annotating Reference (top) Sequence...'

print 'Adding repeats...'
###############################################################################
# put reading in RM *out file, merging, etc into function below

repeats = process_repeat_file(topRepeatMaskFile)
print 'did sort and the population'
print 'adding repeats'
n = 0
ypos = ref_line + 0.05
for R in repeats:
    rs = R[0]
    re = R[1]
    rd = R[2]
    rt = R[3]
    over = 1 # has to be on this chrm at this point, so just force it
    # may change later if we want to draw subset of coords
    if over == 1:
        if rt in color:
            c = color[rt]
        else:
            c = color['Other']
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        if rd == '+':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = re - arrowEndBp
                outFile.write('%f %i %i exon\n' % (ypos,rs,ar))
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,ar,re))
        elif rd == '-':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = rs + arrowEndBp
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,ar))
                outFile.write('%f %i %i exon\n' % (ypos,ar,re))
        else:
            print 'repeat error!  what dir???'
            print R
            sys.exit()
outFile.write('Black\n')
tpos = ypos + 0.02

#set font
outFile.write('/Arial findfont 9 scalefont setfont\n')
outFile.write('%f (Repeats) printname_right\n' % (tpos))
###############################################################################
####################################################
gaps = process_gap_file(data['genomeFragGapsFileName'])
print 'adding contig gaps'
n = 0
ypos = ref_line+0.10
for R in gaps:
    rs = R[0]
    re = R[1]
    over = 1 # has to be on this chrm at this point, so just force it
    # may change later if we want to draw subset of coords
    if over == 1:
        c = 'Black'
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        outFile.write('%f %i %i exon\n' % (ypos,rs,re))
        n += 1
outFile.write('Black\n')
tpos = ypos + 0.02
outFile.write('%f (Gaps) printname_right\n' % (tpos))

####################################################



# would do DUPMASK here, if we had it

###############################################################################
# print out the repeats key, on the left
m = ref_line + 2*0.05 #have two rows...
ys = m + 0.1 # print key 2 lines above htis
outFile.write('/Arial-Bold findfont 10 scalefont setfont\n')
outFile.write('/yh drop %f mul def\n' % (ys)) # transform to right coordinates
# print key at yh
outFile.write('0 yh moveto\n')
outFile.write('Black\n')
outFile.write('(Other     ) 100 string cvs show\n')
outFile.write('DkGray\n')
outFile.write('(Simple Repeat     ) 100 string cvs show\n')
outFile.write('LtGray\n\n')
outFile.write('(Low Complexity     ) 100 string cvs show\n')
outFile.write('Pink\n')
outFile.write('(DNA     ) 100 string cvs show\n')
outFile.write('Orange\n')
outFile.write('(LTR     ) 100 string cvs show\n')
outFile.write('Green\n')
outFile.write('(LINE     ) 100 string cvs show\n')
outFile.write('Purple\n')
outFile.write('(SINE     ) 100 string cvs show\n')
outFile.write('Black\n')
###############################################################################

#Now annotating the fosmid allele
print 'Now annotating the clone (bottom)'
repeats = process_repeat_file(bottomRepeatMaskFile)
print 'Did sort'
print 'Adding repeats....'
print 'adding repeats'
n = 0
ypos = fos_line-0.10
for R in repeats:
    rs = R[0]
    re = R[1]
    rd = R[2]
    rt = R[3]
    over = 1 # has to be on this chrm at this point, so just force it
    # may change later if we want to draw subset of coords
    if over == 1:
        if rt in color:
            c = color[rt]
        else:
            c = color['Other']
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        if rd == '+':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = re - arrowEndBp
                outFile.write('%f %i %i exon\n' % (ypos,rs,ar))
                outFile.write('%f %i %i RIGHT_end\n' % (ypos,ar,re))
        elif rd == '-':
            if eLen < arrowEndBp:
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,re))
                n += 1
            else:
                ar = rs + arrowEndBp
                outFile.write('%f %i %i LEFT_end\n' % (ypos,rs,ar))
                outFile.write('%f %i %i exon\n' % (ypos,ar,re))
        else:
            print 'repeat error!  what dir???'
            print R
            sys.exit()
outFile.write('Black\n')
tpos = ypos + 0.02
outFile.write('/Arial findfont 9 scalefont setfont\n')
outFile.write('%f (Repeats) printname_right\n' % (tpos))

####################################################
gaps = process_gap_file(data['contigSeqGapsFileName'])
print 'adding contig gaps'
n = 0
ypos = fos_line-0.15
for R in gaps:
    rs = R[0]
    re = R[1]
    over = 1 # has to be on this chrm at this point, so just force it
    # may change later if we want to draw subset of coords
    if over == 1:
        c = 'Black'
        eLen = re - rs  + 1
        outFile.write('%s\n' % c)
        outFile.write('%f %i %i exon\n' % (ypos,rs,re))
        n += 1
outFile.write('Black\n')
tpos = ypos + 0.02
outFile.write('%f (Gaps) printname_right\n' % (tpos))

####################################################


# Dupmask would go here if we were doing it for the fosmid


######ADD IN THE BREAK POINTS#####################
# print in the breakpoints here now....

#draw the left bp

genomeP = data['leftBpChromCoords'] - data['chromFragStart'] + 1
contigP = data['leftBpContigCoord']

outFile.write('Blue\n')
outFile.write('%i %f %i %f line_draw_thin\n' % (contigP,fos_line,genomeP,ref_line))

genomeP = data['rightBpChromCoords'] - data['chromFragStart'] + 1
contigP = data['rightBpContigCoord']
outFile.write('Red\n')
outFile.write('%i %f %i %f line_draw_thin\n' % (contigP,fos_line,genomeP,ref_line))
outFile.write('Black\n')


y = fos_line + 0.02
cs = data['contigLeftFragAlignBegin']
ce = data['contigLeftFragAlignEnd']
outFile.write('Blue\n')
outFile.write('%i %f %i %f line_draw\n' % (cs,y,ce,y))

y = fos_line + 0.03
cs = data['contigRightFragAlignBegin']
ce = data['contigRightFragAlignEnd']
outFile.write('Red\n')
outFile.write('%i %f %i %f line_draw\n' % (cs,y,ce,y))

y = ref_line - 0.03
cs = data['genomeFragAlignBegin'] - data['chromFragStart'] + 1
ce = data['genomeFragAlignEnd'] - data['chromFragStart'] + 1
outFile.write('Black\n')
outFile.write('%i %f %i %f line_draw\n' % (cs,y,ce,y))



#    leftSeq = contigSeq[data['contigLeftFragAlignBegin']-1:data['contigLeftFragAlignEnd']]
#    rightSeq = contigSeq[data['contigRightFragAlignBegin']-1:data['contigRightFragAlignEnd']]
#    region = data['chromName'] + ':' + str(data['genomeFragAlignBegin']) + '-' + str(data['genomeFragAlignEnd'])


# print out the sequences used in the alignment


######

#******************************************************
outFile.write(line) # last line we didn't do yet, should be blank
outFile.write('Black\n')

# set scale bar at the appropiate position

outFile.write('/y1 drop %f mul def\n' % fos_line)
outFile.write('/ybot y1 -0.3 cm add def\n')
outFile.write('/ybot2 y1 -0.5 cm add def\n')

while True:
    l = inFile.readline()
    if l == '':
        break
    # replace hard coded positions with pos relative to the fos_line
    l = l.replace('-0.3 cm','ybot')
    l = l.replace('-.3 cm','ybot')

    l = l.replace('-0.5 cm','ybot2')
    l = l.replace('-.5 cm','ybot2')
    outFile.write(l)


inFile.close()
outFile.close()

print 'DONE with annotation adding!'
# ok, now I think we may be done....

