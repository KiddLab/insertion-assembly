import genutils
import os.path
from optparse import  OptionParser
#############################################################################

USAGE = """python combine-pdf-images.py   
                                  --in <file of 3way align PDFs to merge>
                                  --output <file name of merged PDFs>
                                  --miropeats_base <miropeats base dir>
                                  
        combined PDF images into single PDF for review


"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inputFileName', help = 'file of initial miropeats input')
parser.add_option('--out',dest='outFileName', help = 'file of merged PDFs ')
parser.add_option('--miropeats_base',dest='mrDir', help = 'directory of miropeats files')


(options,args)=parser.parse_args()

if options.inputFileName is None:
    parser.error('input miropeats result file not given')
if options.mrDir is None:
    parser.error('miropeats dir not given')
if options.outFileName is None:
    parser.error('output merged files not given')

if options.outFileName[-4:] != '.pdf':
    parser.error('output file name should end in .pdf')
        
#############################################################################
if options.mrDir[-1] != '/':
     options.mrDir  += '/'

baseDir = options.mrDir


have = 0
noHave = 0

toCombine = []
inFileName = options.inputFileName
inFile = open(inFileName,'r')
for line in inFile:
    line = line.rstrip()
    line = line.split()
    siteID = line[0]
    c = line[1]
    if siteID == 'siteID':
        continue
    
    toCombine.append([siteID,c])
inFile.close()

# 'change _ to '.' in names to sort

print 'Have %i file names ' % (len(toCombine))
for i in range(len(toCombine)):
    toCombine[i][0] = toCombine[i][0].replace('_','.')

toCombine.sort()

for i in range(len(toCombine)):
    toCombine[i][0] = toCombine[i][0].replace('.','_')
    s = baseDir + toCombine[i][1] + '/' +  toCombine[i][0] + '/' +  toCombine[i][0] + '.3way.combined.annotated.pdf'
    toCombine[i] = s


newFileName = options.outFileName

tmpFileName = newFileName + '.tmp.pdf'

# prime the pump
cmd = 'cp %s %s' % (toCombine[0],tmpFileName)
print cmd
genutils.runCMD(cmd)
cmd = 'cp %s %s' % (tmpFileName,newFileName)
print cmd
genutils.runCMD(cmd)


for i in range(1,len(toCombine)):
    cmd = 'cpdf %s %s -o %s' % (tmpFileName,toCombine[i],newFileName)
    print cmd
    genutils.runCMD(cmd)
    cmd = 'cp %s %s' % (newFileName,tmpFileName)
    genutils.runCMD(cmd)

cmd = 'rm %s' % tmpFileName
print cmd
genutils.runCMD(cmd)