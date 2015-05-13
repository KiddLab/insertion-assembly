import sys
from optparse import  OptionParser

##########################################################################################
def check_hasMin(hits):
    for i in hits:
        if i[3] >= minRMsize:
            return True
    return False
##########################################################################################
def check_hasMin_OneSide(hits):
    for i in hits:
        if i[3] >= minRMsize and ((i[1] >= minSide) or (i[2] >= minSide)):
            return True
    return False
##########################################################################################
def check_hasMin_TwoSides(hits):
    for i in hits:
        if i[3] >= minRMsize and (i[1] >= minSide) and (i[2] >= minSide):
            return True
    return False
##########################################################################################
#############################################################################

USAGE = """09_select-longest.py --in <selected full fragment>

"""
parser = OptionParser(USAGE)
parser.add_option('--in',dest='inFileName', help = 'sel full fragment')


(options,args)=parser.parse_args()

if options.inFileName is None:
    parser.error('file name for parsed assembly table not given')
    
#############################################################################



minRMsize = 30
minSide = 30

outFileName = options.inFileName + '.longest'
print 'writing output to',outFileName


inFile = open(options.inFileName,'r')
outFile = open(outFileName,'w')

for line in inFile:
    line = line.rstrip()
    line = line.split()
    if line[1] == '1':
        nl = '\t'.join(line) + '\n'
        outFile.write(nl)
        continue
    c_lens = []
    contigs = line[2:]
    for con in contigs:
        t = con.split(':')
        l = int(t[1]) + int(t[2])
        c_lens.append(l)
    m = max(c_lens)
    nl = [line[0],line[1]]
    for i in range(len(c_lens)):
        if c_lens[i] == m:
            nl.append(contigs[i])
            break
    nl = '\t'.join(nl) + '\n'
    outFile.write(nl)

outFile.close()
inFile.close()